#!/usr/bin/env python3
"""
erbb2_sv_compare.py  — SvABA DP 필드 해석 수정판

SvABA sv.vcf FORMAT 필드 해석:
  AD = alt-supporting reads (SR + DR)
  DP = ref-supporting reads   ← total depth가 *아님*!
  total depth = AD + DP
  VAF = AD / (AD + DP)

이 수정 전에는 VAF = AD / DP 로 계산해서 AD > DP인 경우
VAF > 1 이 되는 버그가 있었음.
"""

import argparse
import csv
import gzip
import re
import sys
from collections import OrderedDict

BND_ALT_RE = re.compile(r'[\[\]]([^\[\]:]+):(\d+)[\[\]]')


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_region(region_str):
    chrom, span = region_str.split(":")
    start, end = span.replace(",", "").split("-")
    return chrom, int(start), int(end)


def parse_vcf_breakpoints(path):
    """
    SvABA sv.vcf를 읽어서 breakpoint pair 단위로 정리.

    SvABA FORMAT 해석:
      AD = alt reads (SR + DR)
      DP = ref reads  (≠ total depth)
      total_depth = AD + DP
      VAF = AD / total_depth
    """
    records = {}

    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
            fmt  = fields[8] if len(fields) > 8 else ""
            sample = fields[9] if len(fields) > 9 else ""

            info_map = {}
            for kv in info.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info_map[k] = v
                else:
                    info_map[kv] = True

            mate_match = BND_ALT_RE.search(alt)
            if not mate_match:
                continue
            mate_chrom, mate_pos = mate_match.group(1), int(mate_match.group(2))

            fmt_keys   = fmt.split(":")
            sample_vals = sample.split(":")
            sample_map = dict(zip(fmt_keys, sample_vals))

            def to_int(x):
                try:
                    return int(x)
                except (ValueError, TypeError):
                    return 0

            ad = to_int(sample_map.get("AD", 0))
            dp_ref = to_int(sample_map.get("DP", 0))   # ref reads (SvABA 정의)
            sr = to_int(sample_map.get("SR", 0))
            dr = to_int(sample_map.get("DR", 0))

            total_depth = ad + dp_ref                   # 실제 total depth
            vaf = (ad / total_depth) if total_depth > 0 else float("nan")

            records[vid] = {
                "id": vid,
                "chrom": chrom,
                "pos": int(pos),
                "mate_chrom": mate_chrom,
                "mate_pos": mate_pos,
                "mateid": info_map.get("MATEID"),
                "ad": ad,
                "dp": total_depth,   # 이후 코드에서는 total depth로 통일
                "sr": sr,
                "dr": dr,
                "vaf": vaf,
            }

    # mate pair → 하나의 junction
    seen = set()
    junctions = []
    for vid, rec in records.items():
        if vid in seen:
            continue
        mateid = rec["mateid"]
        seen.add(vid)
        if mateid and mateid in records:
            seen.add(mateid)

        a = (rec["chrom"], rec["pos"])
        b = (rec["mate_chrom"], rec["mate_pos"])
        if a <= b:
            chrom1, pos1, chrom2, pos2 = rec["chrom"], rec["pos"], rec["mate_chrom"], rec["mate_pos"]
        else:
            chrom1, pos1, chrom2, pos2 = rec["mate_chrom"], rec["mate_pos"], rec["chrom"], rec["pos"]

        junctions.append({
            "id": vid,
            "chrom1": chrom1, "pos1": pos1,
            "chrom2": chrom2, "pos2": pos2,
            "ad": rec["ad"],
            "dp": rec["dp"],   # total depth
            "sr": rec["sr"],
            "dr": rec["dr"],
            "vaf": rec["vaf"],
        })

    return junctions


def overlaps_region(chrom, pos, region_chrom, start, end):
    return chrom == region_chrom and start <= pos <= end


def filter_erbb2_linked(junctions, region_chrom, start, end):
    out = []
    for j in junctions:
        hit1 = overlaps_region(j["chrom1"], j["pos1"], region_chrom, start, end)
        hit2 = overlaps_region(j["chrom2"], j["pos2"], region_chrom, start, end)
        if hit1 or hit2:
            j2 = dict(j)
            j2["interchromosomal"] = (j["chrom1"] != j["chrom2"])
            j2["erbb2_side"] = "1" if hit1 else "2"
            out.append(j2)
    return out


def find_match(query, candidates, tol):
    qc1, qp1, qc2, qp2 = query["chrom1"], query["pos1"], query["chrom2"], query["pos2"]
    for c in candidates:
        cc1, cp1, cc2, cp2 = c["chrom1"], c["pos1"], c["chrom2"], c["pos2"]
        same_order = (qc1 == cc1 and abs(qp1 - cp1) <= tol and
                      qc2 == cc2 and abs(qp2 - cp2) <= tol)
        swapped    = (qc1 == cc2 and abs(qp1 - cp2) <= tol and
                      qc2 == cc1 and abs(qp2 - cp1) <= tol)
        if same_order or swapped:
            return c
    return None


def load_manifest(path):
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rows.append(r)
    return rows


def resolve_manifest_samples(manifest_path, source_barcode, parental_sample_id,
                              vcf_template, unfiltered_template):
    rows = load_manifest(manifest_path)
    matched = [r for r in rows if r.get("source_barcode") == source_barcode]
    if not matched:
        sys.exit(f"[error] manifest에 source_barcode={source_barcode} 행이 없습니다.")

    if not unfiltered_template:
        unfiltered_template = vcf_template.replace(".svaba.sv.vcf", ".svaba.unfiltered.sv.vcf")

    parental_path = None
    clones = OrderedDict()
    for r in matched:
        aa_barcode = r["aa_barcode"]
        sample_id  = r["sample_id"]
        filt_path  = vcf_template.format(aa_barcode=aa_barcode)
        unfilt_path = unfiltered_template.format(aa_barcode=aa_barcode)

        if sample_id == parental_sample_id:
            if parental_path is not None:
                print(f"[warn] parental 중복 매칭, 마지막({aa_barcode}) 사용", file=sys.stderr)
            parental_path = filt_path
        else:
            clones[sample_id] = (filt_path, unfilt_path)

    if parental_path is None:
        sys.exit(f"[error] sample_id={parental_sample_id} (parental)을 못 찾았습니다.")

    print(f"[info] manifest: parental=1, clones={len(clones)}", file=sys.stderr)
    return parental_path, clones


def diagnose_zero_match(junctions, region_chrom, start, end):
    chrom_counts = {}
    for j in junctions:
        chrom_counts[j["chrom1"]] = chrom_counts.get(j["chrom1"], 0) + 1
        chrom_counts[j["chrom2"]] = chrom_counts.get(j["chrom2"], 0) + 1

    print("[diag] ERBB2-linked breakpoint 0개. 아래 확인하세요.", file=sys.stderr)
    print(f"[diag] --region chrom: '{region_chrom}'", file=sys.stderr)
    print(f"[diag] VCF chrom 상위 10개: "
          f"{sorted(chrom_counts.items(), key=lambda x: -x[1])[:10]}", file=sys.stderr)

    if region_chrom not in chrom_counts:
        alt = region_chrom[3:] if region_chrom.startswith("chr") else f"chr{region_chrom}"
        if alt in chrom_counts:
            print(f"[diag] !! VCF는 '{alt}' 표기를 씀. --region을 '{alt}:...' 로 바꾸세요.",
                  file=sys.stderr)
        else:
            print(f"[diag] !! VCF에 '{region_chrom}' 자체가 없음. chrom 표기/build 확인요.",
                  file=sys.stderr)
        return

    on_chrom = [j for j in junctions if j["chrom1"] == region_chrom or j["chrom2"] == region_chrom]
    if on_chrom:
        positions = []
        for j in on_chrom:
            if j["chrom1"] == region_chrom: positions.append(j["pos1"])
            if j["chrom2"] == region_chrom: positions.append(j["pos2"])
        print(f"[diag] '{region_chrom}' breakpoint {len(on_chrom)}개 "
              f"({min(positions):,}~{max(positions):,}). "
              f"--region({start:,}-{end:,})과 좌표 불일치. genome build 확인.",
              file=sys.stderr)


def fmt_vaf(vaf):
    import math
    if isinstance(vaf, float) and not math.isnan(vaf):
        return f"{vaf:.3f}"
    return "NA"


def build_long_rows(erbb2_junctions, clone_data, tol):
    rows = []
    for j in erbb2_junctions:
        rows.append({
            "bp_id": j["id"], "chrom1": j["chrom1"], "pos1": j["pos1"],
            "chrom2": j["chrom2"], "pos2": j["pos2"],
            "interchromosomal": j["interchromosomal"], "erbb2_side": j["erbb2_side"],
            "sample_id": "parental", "sample_type": "parental",
            "status": "source", "AD": j["ad"], "DP": j["dp"],
            "VAF": fmt_vaf(j["vaf"]),
        })
        for name, (filt_j, unfilt_j) in clone_data.items():
            m = find_match(j, filt_j, tol)
            status = "detected_filtered"
            if m is None:
                m = find_match(j, unfilt_j, tol)
                status = "detected_unfiltered" if m is not None else "not_detected"

            if m is not None:
                ad, dp = m["ad"], m["dp"]
                vaf = fmt_vaf(m["vaf"])
            else:
                ad, dp, vaf = "NA", "NA", "NA"

            rows.append({
                "bp_id": j["id"], "chrom1": j["chrom1"], "pos1": j["pos1"],
                "chrom2": j["chrom2"], "pos2": j["pos2"],
                "interchromosomal": j["interchromosomal"], "erbb2_side": j["erbb2_side"],
                "sample_id": name, "sample_type": "clone",
                "status": status, "AD": ad, "DP": dp, "VAF": vaf,
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parental")
    ap.add_argument("--clone", action="append", default=[])

    ap.add_argument("--manifest")
    ap.add_argument("--source-barcode")
    ap.add_argument("--parental-sample-id", default="parental")
    ap.add_argument("--vcf-template")
    ap.add_argument("--unfiltered-template")

    ap.add_argument("--region", default="chr17:39687914-39730426")
    ap.add_argument("--flank", type=int, default=2_000_000)
    ap.add_argument("--tol", type=int, default=500)
    ap.add_argument("--format", choices=["long", "wide", "both"], default="long")
    ap.add_argument("--out", default="erbb2_sv_long.tsv")
    args = ap.parse_args()

    region_chrom, start, end = parse_region(args.region)
    start = max(0, start - args.flank)
    end   = end + args.flank

    if args.manifest:
        if not args.source_barcode or not args.vcf_template:
            sys.exit("[error] --manifest 모드에서는 --source-barcode 와 --vcf-template 필수.")
        parental_path, clones = resolve_manifest_samples(
            args.manifest, args.source_barcode, args.parental_sample_id,
            args.vcf_template, args.unfiltered_template,
        )
    else:
        if not args.parental or not args.clone:
            sys.exit("[error] --parental 과 --clone(1개 이상) 필요.")
        parental_path = args.parental
        clones = OrderedDict()
        for c in args.clone:
            name, paths = c.split("=", 1)
            if ":" in paths:
                filt_path, unfilt_path = paths.split(":", 1)
            else:
                filt_path, unfilt_path = paths, None
            clones[name] = (filt_path, unfilt_path)

    print(f"[info] parental: {parental_path}", file=sys.stderr)
    parental_junctions = parse_vcf_breakpoints(parental_path)
    erbb2_junctions    = filter_erbb2_linked(parental_junctions, region_chrom, start, end)
    print(f"[info] parental breakpoints total={len(parental_junctions)}, "
          f"ERBB2-linked(±{args.flank}bp)={len(erbb2_junctions)}", file=sys.stderr)

    if len(erbb2_junctions) == 0:
        diagnose_zero_match(parental_junctions, region_chrom, start, end)

    clone_data = {}
    for name, (filt_path, unfilt_path) in clones.items():
        filt_j   = parse_vcf_breakpoints(filt_path)
        unfilt_j = parse_vcf_breakpoints(unfilt_path) if unfilt_path else []
        clone_data[name] = (filt_j, unfilt_j)
        print(f"[info] clone {name}: filtered={len(filt_j)} unfiltered={len(unfilt_j)}",
              file=sys.stderr)

    long_header = ["bp_id", "chrom1", "pos1", "chrom2", "pos2", "interchromosomal",
                   "erbb2_side", "sample_id", "sample_type", "status", "AD", "DP", "VAF"]
    long_rows = build_long_rows(erbb2_junctions, clone_data, args.tol)

    wide_header = ["bp_id", "chrom1", "pos1", "chrom2", "pos2",
                   "interchromosomal", "erbb2_side",
                   "parental_AD", "parental_DP", "parental_VAF"]
    for name in clones:
        wide_header += [f"{name}_status", f"{name}_AD", f"{name}_DP", f"{name}_VAF"]

    wide_rows = []
    for j in erbb2_junctions:
        row = [j["id"], j["chrom1"], j["pos1"], j["chrom2"], j["pos2"],
               j["interchromosomal"], j["erbb2_side"],
               j["ad"], j["dp"], fmt_vaf(j["vaf"])]
        for name, (filt_j, unfilt_j) in clone_data.items():
            m = find_match(j, filt_j, args.tol)
            status = "detected_filtered"
            if m is None:
                m = find_match(j, unfilt_j, args.tol)
                status = "detected_unfiltered" if m is not None else "not_detected"
            if m is not None:
                row += [status, m["ad"], m["dp"], fmt_vaf(m["vaf"])]
            else:
                row += [status, "NA", "NA", "NA"]
        wide_rows.append(row)

    def write_long(path):
        with open(path, "w") as fh:
            fh.write("\t".join(long_header) + "\n")
            for r in long_rows:
                fh.write("\t".join(str(r[c]) for c in long_header) + "\n")

    def write_wide(path):
        with open(path, "w") as fh:
            fh.write("\t".join(wide_header) + "\n")
            for r in wide_rows:
                fh.write("\t".join(str(x) for x in r) + "\n")

    if args.format == "long":
        write_long(args.out)
        print(f"[done] long, {len(long_rows)} rows -> {args.out}", file=sys.stderr)
    elif args.format == "wide":
        write_wide(args.out)
        print(f"[done] wide, {len(wide_rows)} rows -> {args.out}", file=sys.stderr)
    else:
        long_path = args.out
        wide_path = (args.out.rsplit(".", 1)[0] + ".wide.tsv"
                     if "." in args.out else args.out + ".wide")
        write_long(long_path)
        write_wide(wide_path)
        print(f"[done] long -> {long_path} ({len(long_rows)} rows), "
              f"wide -> {wide_path} ({len(wide_rows)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
