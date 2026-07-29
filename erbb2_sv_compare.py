#!/usr/bin/env python3
"""
erbb2_sv_compare.py

Parental SvABA sv.vcf 기준으로 ERBB2-linked breakpoint를 뽑고,
각 clone의 sv.vcf(+필요시 unfiltered.sv.vcf)에서 같은 breakpoint가
검출되는지/AD,DP,VAF가 어떻게 다른지 매트릭스로 만든다.

사용 예 1 (수동으로 clone 하나씩 지정, 클론 적을 때):
python erbb2_sv_compare.py \
  --parental /path/parental/ECTRES-H2170-...-parental.svaba.sv.vcf \
  --clone CLONE1=/path/clone1/ECTRES-...svaba.sv.vcf:/path/clone1/ECTRES-...svaba.unfiltered.sv.vcf \
  --clone CLONE2=/path/clone2/ECTRES-...svaba.sv.vcf:/path/clone2/ECTRES-...svaba.unfiltered.sv.vcf \
  --region chr17:39687914-39730426 \
  --flank 2000000 \
  --tol 500 \
  --out erbb2_sv_matrix.tsv

사용 예 2 (manifest로 클론 많을 때 자동 처리):
sample_map.csv 형식: aa_barcode,source_barcode,sample_id
  ECTRES-H2170-0001-TPX-A09-WGS-6XP120,H2170,NCI_9
  ECTRES-H2170-0001-TPX-A01-WGS-XXXXXX,H2170,parental
  ...

python erbb2_sv_compare.py \
  --manifest sample_map.csv \
  --source-barcode H2170 \
  --parental-sample-id parental \
  --vcf-template "/mnt/NAS3/home/jiwon/HL-NF/scratch/ECTRES/results/svaba/svaba_ss/{aa_barcode}/vcfs/{aa_barcode}.svaba.sv.vcf" \
  --region chr17:39687914-39730426 \
  --flank 2000000 \
  --tol 500 \
  --format long \
  --out erbb2_sv_long.tsv

  - {aa_barcode}는 manifest의 aa_barcode 컬럼 값으로 자동 치환됨
  - unfiltered vcf 경로는 --unfiltered-template을 안 주면
    --vcf-template에서 ".svaba.sv.vcf" -> ".svaba.unfiltered.sv.vcf" 로 자동 변환해서 씀
    (실제 폴더 구조가 그 명명규칙을 따른다는 전제, 다르면 --unfiltered-template으로 직접 지정)
  - manifest에서 source_barcode == --source-barcode 인 행 중
    sample_id == --parental-sample-id (기본값 "parental") 인 행을 parental로,
    나머지 전부를 clone으로 자동 등록함

출력 포맷:
  --format long (기본값, 권장): 한 행 = breakpoint x sample 하나.
    컬럼: bp_id, chrom1, pos1, chrom2, pos2, interchromosomal, erbb2_side,
          sample_id, sample_type(parental/clone), status, AD, DP, VAF
    -> R/Python에서 plotting(heatmap, dotplot 등) 하기 편한 tidy format.
       클론이 몇 개든 컬럼 수가 안 늘어남.
  --format wide: 기존 방식. clone마다 4개 컬럼(_status/_AD/_DP/_VAF)이 옆으로 늘어남.
    클론 적을 때 눈으로 훑어보기엔 편하지만, 많아지면 컬럼이 폭발함.
  --format both: 두 파일 다 저장 (long은 --out 그대로, wide는 --out에 .wide.tsv 붙여서 저장)

주의:
- ERBB2 좌표는 기본값을 GRCh38로 넣어뒀음. hg19면 --region으로 직접 덮어쓸 것
  (hg19 ERBB2 ~ chr17:37,844,167-37,886,679)
- --flank 는 ERBB2 유전자 자체가 아니라 amplicon 전체를 보기 위한 여유 범위.
  AA로 정의된 amplicon segment bed가 있으면 그 범위로 더 정확히 좁힐 수 있음
  (--region을 amplicon min-max 좌표로 바꿔서 사용).
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
    각 BND record는 mate record와 짝을 이루므로, MATEID로 짝지어
    하나의 junction당 하나의 row만 반환 (canonical: (chromA,posA) <= (chromB,posB)).

    return: list of dict with keys:
      id, chrom1, pos1, chrom2, pos2, ad, dp, sr, dr, vaf
    """
    records = {}  # id -> dict(chrom,pos,alt,mateid,ad,dp,sr,dr)
    format_idx = None

    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
            fmt = fields[8] if len(fields) > 8 else ""
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
                # BND가 아니면 skip (svaba sv.vcf는 전부 BND여야 정상)
                continue
            mate_chrom, mate_pos = mate_match.group(1), int(mate_match.group(2))

            fmt_keys = fmt.split(":")
            sample_vals = sample.split(":")
            sample_map = dict(zip(fmt_keys, sample_vals))

            def to_int(x):
                try:
                    return int(x)
                except (ValueError, TypeError):
                    return 0

            ad = to_int(sample_map.get("AD", 0))
            dp = to_int(sample_map.get("DP", 0))
            sr = to_int(sample_map.get("SR", 0))
            dr = to_int(sample_map.get("DR", 0))

            records[vid] = {
                "id": vid,
                "chrom": chrom,
                "pos": int(pos),
                "mate_chrom": mate_chrom,
                "mate_pos": mate_pos,
                "mateid": info_map.get("MATEID"),
                "ad": ad,
                "dp": dp,
                "sr": sr,
                "dr": dr,
            }

    # mate pair를 하나의 junction으로 합치기
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

        dp = rec["dp"]
        ad = rec["ad"]
        vaf = (ad / dp) if dp > 0 else float("nan")

        junctions.append({
            "id": vid,
            "chrom1": chrom1, "pos1": pos1,
            "chrom2": chrom2, "pos2": pos2,
            "ad": ad, "dp": dp, "sr": rec["sr"], "dr": rec["dr"],
            "vaf": vaf,
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
    """query breakpoint와 candidates(list of junction dict) 중
    좌표가 tol 이내로 맞는 (순서 무관) 첫 매치를 반환, 없으면 None."""
    qc1, qp1, qc2, qp2 = query["chrom1"], query["pos1"], query["chrom2"], query["pos2"]
    for c in candidates:
        cc1, cp1, cc2, cp2 = c["chrom1"], c["pos1"], c["chrom2"], c["pos2"]
        same_order = (qc1 == cc1 and abs(qp1 - cp1) <= tol and
                      qc2 == cc2 and abs(qp2 - cp2) <= tol)
        swapped = (qc1 == cc2 and abs(qp1 - cp2) <= tol and
                   qc2 == cc1 and abs(qp2 - cp1) <= tol)
        if same_order or swapped:
            return c
    return None


def load_manifest(path):
    """sample_map.csv (aa_barcode,source_barcode,sample_id) 읽기"""
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for r in reader:
            rows.append(r)
    return rows


def resolve_manifest_samples(manifest_path, source_barcode, parental_sample_id,
                              vcf_template, unfiltered_template):
    """manifest에서 source_barcode가 일치하는 행들을 골라
    parental 경로 1개 + clone 경로 dict로 분리해서 반환."""
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
        sample_id = r["sample_id"]
        filt_path = vcf_template.format(aa_barcode=aa_barcode)
        unfilt_path = unfiltered_template.format(aa_barcode=aa_barcode)

        if sample_id == parental_sample_id:
            if parental_path is not None:
                print(f"[warn] parental sample_id={parental_sample_id} 가 두 번 이상 매칭됨, "
                      f"마지막 것({aa_barcode})을 사용", file=sys.stderr)
            parental_path = filt_path
        else:
            clones[sample_id] = (filt_path, unfilt_path)

    if parental_path is None:
        sys.exit(f"[error] manifest에서 sample_id={parental_sample_id} (parental)을 못 찾았습니다. "
                  f"--parental-sample-id 로 정확한 값을 지정하세요.")

    print(f"[info] manifest에서 source_barcode={source_barcode}: "
          f"parental=1, clones={len(clones)}", file=sys.stderr)
    return parental_path, clones


def diagnose_zero_match(junctions, region_chrom, start, end):
    """ERBB2-linked breakpoint가 0개일 때, chrom naming/build 문제인지
    힌트를 주기 위한 진단 정보 출력."""
    chrom_counts = {}
    for j in junctions:
        chrom_counts[j["chrom1"]] = chrom_counts.get(j["chrom1"], 0) + 1
        chrom_counts[j["chrom2"]] = chrom_counts.get(j["chrom2"], 0) + 1

    print("[diag] ERBB2-linked breakpoint가 0개입니다. 아래를 확인하세요.", file=sys.stderr)
    print(f"[diag] --region 에 쓴 chrom 이름: '{region_chrom}'", file=sys.stderr)
    print(f"[diag] VCF에 실제 등장하는 chrom 이름들(상위 10개): "
          f"{sorted(chrom_counts.items(), key=lambda x: -x[1])[:10]}", file=sys.stderr)

    if region_chrom not in chrom_counts:
        alt = region_chrom[3:] if region_chrom.startswith("chr") else f"chr{region_chrom}"
        if alt in chrom_counts:
            print(f"[diag] !! VCF는 '{alt}' 표기를 쓰는데 --region은 '{region_chrom}'을 씀. "
                  f"--region 값을 '{alt}:...' 형태로 바꿔서 재시도하세요.", file=sys.stderr)
        else:
            print(f"[diag] !! VCF에 '{region_chrom}' 자체가 한 번도 안 나옴. "
                  f"chrom 이름 표기나 reference build를 다시 확인하세요.", file=sys.stderr)
        return

    on_chrom = [j for j in junctions if j["chrom1"] == region_chrom or j["chrom2"] == region_chrom]
    if on_chrom:
        positions = []
        for j in on_chrom:
            if j["chrom1"] == region_chrom:
                positions.append(j["pos1"])
            if j["chrom2"] == region_chrom:
                positions.append(j["pos2"])
        print(f"[diag] '{region_chrom}' 위에는 breakpoint {len(on_chrom)}개 있음 "
              f"(좌표 범위: {min(positions):,} ~ {max(positions):,}). "
              f"--region({start:,}-{end:,})과 안 겹치는 것 같으니, "
              f"genome build(hg19 vs GRCh38)가 다를 가능성이 높습니다. "
              f"--flank를 늘려보거나 --region 좌표를 직접 다시 확인하세요.", file=sys.stderr)


def build_long_rows(erbb2_junctions, clone_data, tol):
    """long(tidy) format row 생성: 한 행 = breakpoint x sample"""
    rows = []
    for j in erbb2_junctions:
        # parental 자기 자신도 한 row로 포함 (sample_type=parental)
        rows.append({
            "bp_id": j["id"], "chrom1": j["chrom1"], "pos1": j["pos1"],
            "chrom2": j["chrom2"], "pos2": j["pos2"],
            "interchromosomal": j["interchromosomal"], "erbb2_side": j["erbb2_side"],
            "sample_id": "parental", "sample_type": "parental",
            "status": "source", "AD": j["ad"], "DP": j["dp"],
            "VAF": f"{j['vaf']:.3f}" if j["dp"] > 0 else "NA",
        })
        for name, (filt_j, unfilt_j) in clone_data.items():
            m = find_match(j, filt_j, tol)
            status = "detected_filtered"
            if m is None:
                m = find_match(j, unfilt_j, tol)
                status = "detected_unfiltered" if m is not None else "not_detected"

            if m is not None:
                ad, dp = m["ad"], m["dp"]
                vaf = f"{m['vaf']:.3f}" if dp > 0 else "NA"
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
    ap.add_argument("--parental", help="parental svaba.sv.vcf 경로 (manifest 안 쓸 때)")
    ap.add_argument("--clone", action="append", default=[],
                     help="NAME=filtered_vcf[:unfiltered_vcf] 형식, 여러 번 지정 가능 (manifest 안 쓸 때)")

    ap.add_argument("--manifest", help="sample_map.csv 경로 (aa_barcode,source_barcode,sample_id)")
    ap.add_argument("--source-barcode", help="manifest에서 필터할 source_barcode 값 (예: H2170)")
    ap.add_argument("--parental-sample-id", default="parental",
                     help="manifest에서 parental로 취급할 sample_id 값 (기본 'parental')")
    ap.add_argument("--vcf-template",
                     help="manifest 모드일 때 filtered sv.vcf 경로 템플릿. {aa_barcode} 플레이스홀더 사용")
    ap.add_argument("--unfiltered-template",
                     help="manifest 모드일 때 unfiltered sv.vcf 경로 템플릿. "
                          "안 주면 --vcf-template의 '.svaba.sv.vcf'를 '.svaba.unfiltered.sv.vcf'로 자동 치환")

    ap.add_argument("--region", default="chr17:39687914-39730426",
                     help="ERBB2 좌표 (기본 GRCh38). hg19면 직접 지정")
    ap.add_argument("--flank", type=int, default=2_000_000,
                     help="ERBB2 좌우 추가 검색 범위 (amplicon 전체를 보기 위함, default 2Mb)")
    ap.add_argument("--tol", type=int, default=500,
                     help="breakpoint 좌표 매칭 허용 오차 (bp, default 500)")
    ap.add_argument("--format", choices=["long", "wide", "both"], default="long",
                     help="출력 포맷 (기본 long, plotting용으로 권장)")
    ap.add_argument("--out", default="erbb2_sv_long.tsv")
    args = ap.parse_args()

    region_chrom, start, end = parse_region(args.region)
    start = max(0, start - args.flank)
    end = end + args.flank

    if args.manifest:
        if not args.source_barcode or not args.vcf_template:
            sys.exit("[error] --manifest 모드에서는 --source-barcode 와 --vcf-template 가 필수입니다.")
        parental_path, clones = resolve_manifest_samples(
            args.manifest, args.source_barcode, args.parental_sample_id,
            args.vcf_template, args.unfiltered_template,
        )
    else:
        if not args.parental or not args.clone:
            sys.exit("[error] --manifest를 안 쓸 경우 --parental 과 --clone(1개 이상)이 필요합니다.")
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
    erbb2_junctions = filter_erbb2_linked(parental_junctions, region_chrom, start, end)
    print(f"[info] parental breakpoints total={len(parental_junctions)}, "
          f"ERBB2-linked(±{args.flank}bp)={len(erbb2_junctions)}", file=sys.stderr)

    if len(erbb2_junctions) == 0:
        diagnose_zero_match(parental_junctions, region_chrom, start, end)
        # 0개라도 빈 파일은 만들어둠 (디버깅 후 재실행 흐름이 끊기지 않게)

    clone_data = {}
    for name, (filt_path, unfilt_path) in clones.items():
        filt_j = parse_vcf_breakpoints(filt_path)
        unfilt_j = parse_vcf_breakpoints(unfilt_path) if unfilt_path else []
        clone_data[name] = (filt_j, unfilt_j)
        print(f"[info] clone {name}: filtered={len(filt_j)} unfiltered={len(unfilt_j)}",
              file=sys.stderr)

    long_header = ["bp_id", "chrom1", "pos1", "chrom2", "pos2", "interchromosomal",
                   "erbb2_side", "sample_id", "sample_type", "status", "AD", "DP", "VAF"]
    long_rows = build_long_rows(erbb2_junctions, clone_data, args.tol)

    wide_header = ["bp_id", "chrom1", "pos1", "chrom2", "pos2", "interchromosomal", "erbb2_side",
                   "parental_AD", "parental_DP", "parental_VAF"]
    for name in clones:
        wide_header += [f"{name}_status", f"{name}_AD", f"{name}_DP", f"{name}_VAF"]

    wide_rows = []
    for j in erbb2_junctions:
        row = [j["id"], j["chrom1"], j["pos1"], j["chrom2"], j["pos2"],
               j["interchromosomal"], j["erbb2_side"],
               j["ad"], j["dp"], f"{j['vaf']:.3f}" if j["dp"] > 0 else "NA"]
        for name, (filt_j, unfilt_j) in clone_data.items():
            m = find_match(j, filt_j, args.tol)
            status = "detected_filtered"
            if m is None:
                m = find_match(j, unfilt_j, args.tol)
                status = "detected_unfiltered" if m is not None else "not_detected"
            if m is not None:
                vaf_str = f"{m['vaf']:.3f}" if m["dp"] > 0 else "NA"
                row += [status, m["ad"], m["dp"], vaf_str]
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
        print(f"[done] long format, {len(long_rows)} rows -> {args.out}", file=sys.stderr)
    elif args.format == "wide":
        write_wide(args.out)
        print(f"[done] wide format, {len(wide_rows)} rows -> {args.out}", file=sys.stderr)
    else:
        long_path = args.out
        wide_path = args.out.rsplit(".", 1)[0] + ".wide.tsv" if "." in args.out else args.out + ".wide"
        write_long(long_path)
        write_wide(wide_path)
        print(f"[done] long -> {long_path} ({len(long_rows)} rows), "
              f"wide -> {wide_path} ({len(wide_rows)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()