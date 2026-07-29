#!/usr/bin/env python3
"""
annotate_svaba_vcf.py

For ONE sample: annotate its SvABA SV VCF with AA (AmpliconArchitect)
breakpoint-matching info, based on the SCTG (assembled contig) columns in
AA's *_matched_[unfiltered_]svaba_annotated_aa.tsv.

Logic (as specified):
  - For each VCF record's SCTG (INFO field), look it up against every
    contig name appearing anywhere in AA's SCTG_1bp / SCTG_10bp /
    SCTG_100bp / SCTG_500bp columns.
  - If found, tag with the FINEST (smallest) window it was found at:
    1bp > 10bp > 100bp > 500bp (i.e. prefer 1bp over 10bp, etc.)
  - If not found in any of the four columns: AA_SCTG_TIER=NA
  - Also carry through which AA amplicon / cycle / cycle-class that
    contig belongs to (so ecDNA-like ones are directly visible).
  - Write result as a NEW vcf file (original untouched).

This does ONE sample at a time. Combining multiple samples is a separate,
later step -- not handled here.

Usage:
    python3 annotate_svaba_vcf.py \
        --vcf ECTRES-H2170-0001-TPX-A01-WGS-3YV111_svaba_sv.vcf \
        --aa-tsv ECTRES-H2170-0001-TPX-A01-WGS-3YV111_matched_svaba_annotated_aa.tsv \
        --out ECTRES-H2170-0001-TPX-A01-WGS-3YV111_svaba_sv.ecDNA_annotated.vcf
"""
import argparse
import glob
import os
import re
import pandas as pd

TIER_COLS_IN_ORDER = ["SCTG_1bp", "SCTG_10bp", "SCTG_100bp", "SCTG_500bp"]
TIER_LABEL = {"SCTG_1bp": "1bp", "SCTG_10bp": "10bp", "SCTG_100bp": "100bp", "SCTG_500bp": "500bp"}
TIER_RANK = {"SCTG_1bp": 1, "SCTG_10bp": 2, "SCTG_100bp": 3, "SCTG_500bp": 4}

NEW_INFO_HEADER_LINES = [
    '##INFO=<ID=AA_SCTG_TIER,Number=1,Type=String,Description='
    '"Finest AA (AmpliconArchitect) SCTG match window for this SV\'s contig: '
    '1bp, 10bp, 100bp, 500bp, or NA if not matched at any window">',
    '##INFO=<ID=AA_AMPLICON,Number=.,Type=String,Description='
    '"AA amplicon_id(s) whose cycle breakpoint contig matches this SV (only present if AA_SCTG_TIER != NA)">',
    '##INFO=<ID=AA_CYCLE,Number=.,Type=String,Description='
    '"AA cycle number(s) within the amplicon(s) that this SV\'s contig matches">',
    '##INFO=<ID=AA_CYCLECLASS,Number=.,Type=String,Description='
    '"AA CycleClass(es) (e.g. ecDNA-like, BFB-like, Linear, Invalid) for the matching cycle(s)">',
]


def _amplicon_id_from_filename(fn):
    m = re.search(r"amplicon(\d+)", fn)
    return f"amplicon{m.group(1)}" if m else "NA"


def build_contig_lookup(aa_tsv_path):
    """
    Returns dict: contig_name -> {
        'tier': 'SCTG_1bp' (the finest tier this contig was seen at, across the whole file),
        'amplicons': set(...), 'cycles': set(...), 'cycleclasses': set(...)
    }
    """
    df = pd.read_csv(aa_tsv_path, sep="\t", dtype=str)
    has_cycleclass = "CycleClass" in df.columns

    lookup = {}
    for _, row in df.iterrows():
        amplicon = _amplicon_id_from_filename(row.get("file_name", ""))
        cycle = str(row.get("Cycle", "")).replace("Cycle=", "")
        cycleclass = str(row.get("CycleClass", "")).replace("CycleClass=", "") if has_cycleclass else "NA"

        for col in TIER_COLS_IN_ORDER:
            val = row.get(col)
            if pd.isna(val) or val in ("", "NA"):
                continue
            for contig in str(val).split(","):
                contig = contig.strip()
                if not contig or contig == "NA":
                    continue
                entry = lookup.setdefault(contig, {
                    "tier": col, "amplicons": set(), "cycles": set(), "cycleclasses": set()
                })
                if TIER_RANK[col] < TIER_RANK[entry["tier"]]:
                    entry["tier"] = col
                entry["amplicons"].add(amplicon)
                entry["cycles"].add(cycle)
                entry["cycleclasses"].add(cycleclass)
    return lookup


def get_sctg_from_info(info_str):
    for field in info_str.split(";"):
        if field.startswith("SCTG="):
            return field[len("SCTG="):]
    return None


def annotate_record_info(info_str, lookup):
    sctg_val = get_sctg_from_info(info_str)
    best_tier_rank, best_tier_col = None, None
    amplicons, cycles, cycleclasses = set(), set(), set()

    if sctg_val:
        for contig in sctg_val.split(","):
            contig = contig.strip()
            entry = lookup.get(contig)
            if entry is None:
                continue
            if best_tier_rank is None or TIER_RANK[entry["tier"]] < best_tier_rank:
                best_tier_rank = TIER_RANK[entry["tier"]]
                best_tier_col = entry["tier"]
            amplicons |= entry["amplicons"]
            cycles |= entry["cycles"]
            cycleclasses |= entry["cycleclasses"]

    if best_tier_col is None:
        new_fields = "AA_SCTG_TIER=NA"
    else:
        new_fields = (
            f"AA_SCTG_TIER={TIER_LABEL[best_tier_col]};"
            f"AA_AMPLICON={','.join(sorted(amplicons))};"
            f"AA_CYCLE={','.join(sorted(cycles))};"
            f"AA_CYCLECLASS={','.join(sorted(cycleclasses))}"
        )
    return info_str + ";" + new_fields


def annotate_vcf(vcf_in_path, aa_tsv_path, vcf_out_path):
    lookup = build_contig_lookup(aa_tsv_path)

    with open(vcf_in_path) as fin, open(vcf_out_path, "w") as fout:
        chrom_header_written = False
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue
            if line.startswith("#CHROM"):
                for h in NEW_INFO_HEADER_LINES:
                    fout.write(h + "\n")
                fout.write(line)
                chrom_header_written = True
                continue
            fields = line.rstrip("\n").split("\t")
            info_idx = 7
            fields[info_idx] = annotate_record_info(fields[info_idx], lookup)
            fout.write("\t".join(fields) + "\n")

    n_contigs = len(lookup)
    return n_contigs


def find_one(pattern, exclude_substr=None):
    hits = glob.glob(pattern)
    if exclude_substr:
        hits = [h for h in hits if exclude_substr not in os.path.basename(h)]
    if len(hits) > 1:
        print(f"  [WARN] multiple files match {pattern!r}: {hits} -- using first")
    return hits[0] if hits else None


def discover_sample_files(base, sample):
    """
    Expects the layout:
        {base}/{sample}/align_aa_bp_to_svaba/{sample}_matched_svaba_annotated_aa.tsv
        {base}/{sample}/align_aa_bp_to_svaba/{sample}_matched_unfiltered_svaba_annotated_aa.tsv
        {base}/{sample}/svaba_ss/local_assembly/{sample}.svaba.sv.vcf
        {base}/{sample}/svaba_ss/local_assembly/{sample}.svaba.unfiltered.sv.vcf
    """
    sdir = os.path.join(base, sample)
    aa_dir = os.path.join(sdir, "align_aa_bp_to_svaba")
    vcf_dir = os.path.join(sdir, "svaba_ss", "local_assembly")
    return {
        "aa_annotated_filtered": find_one(
            os.path.join(aa_dir, f"{sample}_matched_svaba_annotated_aa.tsv"), exclude_substr="unfiltered"),
        "aa_annotated_unfiltered": find_one(
            os.path.join(aa_dir, f"{sample}_matched_unfiltered_svaba_annotated_aa.tsv")),
        "vcf_filtered": find_one(
            os.path.join(vcf_dir, f"{sample}.svaba.sv.vcf"), exclude_substr="unfiltered"),
        "vcf_unfiltered": find_one(
            os.path.join(vcf_dir, f"{sample}.svaba.unfiltered.sv.vcf")),
    }


def annotate_sample_by_dir(base, sample, out_dir=None):
    """
    Auto-discover a sample's filtered+unfiltered VCF and matching AA tsv under
    the NAS layout, annotate BOTH, and write into
        {base}/{sample}/ecDNA_annotated_vcf/
    (or `out_dir` if given).
    """
    files = discover_sample_files(base, sample)
    out_dir = out_dir or os.path.join(base, sample, "ecDNA_annotated_vcf")
    os.makedirs(out_dir, exist_ok=True)

    results = {}
    pairs = [
        ("filtered", files["vcf_filtered"], files["aa_annotated_filtered"]),
        ("unfiltered", files["vcf_unfiltered"], files["aa_annotated_unfiltered"]),
    ]
    for vcf_type, vcf_path, aa_path in pairs:
        if vcf_path is None or aa_path is None:
            missing = "vcf" if vcf_path is None else "aa-tsv"
            print(f"[SKIP] {sample} ({vcf_type}): missing {missing}")
            continue
        out_path = os.path.join(out_dir, f"{sample}.svaba.{vcf_type}.sv.ecDNA_annotated.vcf")
        n_contigs = annotate_vcf(vcf_path, aa_path, out_path)
        print(f"[OK] {sample} ({vcf_type}): {n_contigs} AA contigs looked up -> {out_path}")
        results[vcf_type] = out_path
    return results


def annotate_all_samples(base, out_dir_name="ecDNA_annotated_vcf", sample_glob="ECTRES-*"):
    """
    Process every sample directory under `base` (matching sample_glob).
    Never stops on one sample's failure -- logs it and moves on, then
    prints a final OK/SKIPPED/ERROR summary.
    """
    sample_dirs = sorted(d for d in glob.glob(os.path.join(base, sample_glob)) if os.path.isdir(d))
    print(f"Found {len(sample_dirs)} sample dirs under {base}")

    summary = {"OK": [], "PARTIAL": [], "SKIPPED": [], "ERROR": []}
    for sdir in sample_dirs:
        sample = os.path.basename(sdir)
        try:
            results = annotate_sample_by_dir(base, sample)
            if len(results) == 2:
                summary["OK"].append(sample)
            elif len(results) == 1:
                summary["PARTIAL"].append(sample)
            else:
                summary["SKIPPED"].append(sample)
        except Exception as e:
            print(f"[ERROR] {sample}: {e}")
            summary["ERROR"].append(sample)

    print("\n=== Batch summary ===")
    for status, samples in summary.items():
        print(f"{status}: {len(samples)}")
        if status in ("ERROR", "SKIPPED", "PARTIAL"):
            for s in samples:
                print(f"    - {s}")
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", help="Input SvABA sv.vcf (explicit-path mode)")
    ap.add_argument("--aa-tsv", help="AA *_matched_svaba_annotated_aa.tsv (explicit-path mode)")
    ap.add_argument("--out", help="Output annotated VCF path (explicit-path mode)")
    ap.add_argument("--base", help="Base results dir (use alone for ALL samples, or with --sample for one)")
    ap.add_argument("--sample", help="Sample dir name under --base (single-sample auto-discovery mode)")
    ap.add_argument("--out-dir", default=None,
                     help="Output dir for single-sample auto-discovery mode "
                          "(default: {base}/{sample}/ecDNA_annotated_vcf)")
    ap.add_argument("--sample-glob", default="ECTRES-*",
                     help="Glob pattern for sample dir names when processing all samples under --base "
                          "(default: 'ECTRES-*')")
    args = ap.parse_args()

    if args.base and args.sample:
        annotate_sample_by_dir(args.base, args.sample, out_dir=args.out_dir)
    elif args.base:
        annotate_all_samples(args.base, sample_glob=args.sample_glob)
    elif args.vcf and args.aa_tsv and args.out:
        n_contigs = annotate_vcf(args.vcf, args.aa_tsv, args.out)
        print(f"AA contig lookup built from {n_contigs} distinct contigs")
        print(f"Wrote annotated VCF: {args.out}")
    else:
        ap.error("Provide either (--vcf --aa-tsv --out) for a single file, "
                  "(--base --sample) for one sample, or (--base) alone for ALL samples under it.")
