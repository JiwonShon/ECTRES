#!/usr/bin/env python3
"""
combine_filtered_vcf_to_csv.py

Find every sample's FILTERED, ecDNA-annotated SvABA VCF
    {base}/{sample}/ecDNA_annotated_vcf/{sample}.svaba.filtered.sv.ecDNA_annotated.vcf
(as produced by annotate_svaba_vcf.py), parse each into a table with:
    - standard VCF columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER
    - every INFO key as its own column (SVTYPE, SCTG, SPAN, MAPQ, AA_SCTG_TIER,
      AA_AMPLICON, AA_CYCLE, AA_CYCLECLASS, ...)
    - every FORMAT key as its own column (GT, AD, DP, GQ, PL, SR, DR, LR, LO)
    - aliquot_barcode: the sample name parsed from the VCF's sample column
      header (e.g. "ECTRES-H2170-0001-TPX-A31-WGS-VZD4GS.realn.mdup.bqsr.bam"
      -> "ECTRES-H2170-0001-TPX-A31-WGS-VZD4GS")

...then concatenates all samples into ONE combined CSV.

Usage:
    python3 combine_filtered_vcf_to_csv.py \
        --base /mnt/NAS3/.../10X \
        --out  /mnt/NAS3/.../10X/ecDNA_SVABA_filtered_combined.csv
"""
import argparse
import glob
import os
import pandas as pd


def _parse_info(info_str):
    d = {}
    for field in info_str.split(";"):
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
        else:
            d[field] = True
    return d


def parse_annotated_vcf(vcf_path):
    """Returns a DataFrame, one row per VCF record, INFO/FORMAT split into columns."""
    header_cols = None
    n_skip = 0
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                n_skip += 1
                continue
            if line.startswith("#CHROM"):
                header_cols = line.lstrip("#").strip().split("\t")
                break
    if header_cols is None:
        raise ValueError(f"No #CHROM header found in {vcf_path}")

    std_cols = {"CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"}
    sample_col = [c for c in header_cols if c not in std_cols][0]
    aliquot_barcode = sample_col.split(".")[0]

    raw = pd.read_csv(vcf_path, sep="\t", skiprows=n_skip + 1, header=None,
                       names=header_cols, dtype=str)

    rows = []
    for _, r in raw.iterrows():
        row = {
            "aliquot_barcode": aliquot_barcode,
            "CHROM": r["CHROM"], "POS": r["POS"], "ID": r["ID"],
            "REF": r["REF"], "ALT": r["ALT"], "QUAL": r["QUAL"], "FILTER": r["FILTER"],
        }
        row.update(_parse_info(r["INFO"]))
        fmt_keys = r["FORMAT"].split(":")
        fmt_vals = r[sample_col].split(":")
        row.update(dict(zip(fmt_keys, fmt_vals)))
        rows.append(row)

    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="Base results dir containing sample subdirs")
    ap.add_argument("--out", required=True, help="Output combined CSV path")
    ap.add_argument("--pattern", default="*/ecDNA_annotated_vcf/*.svaba.filtered.sv.ecDNA_annotated.vcf",
                     help="Glob pattern (relative to --base) for filtered annotated VCFs")
    args = ap.parse_args()

    vcf_files = sorted(glob.glob(os.path.join(args.base, args.pattern)))
    print(f"Found {len(vcf_files)} filtered annotated VCF files")
    if not vcf_files:
        raise SystemExit("No files matched -- check --base / --pattern, "
                          "and that annotate_svaba_vcf.py has been run.")

    dfs = []
    for f in vcf_files:
        try:
            dfs.append(parse_annotated_vcf(f))
        except Exception as e:
            print(f"[WARN] failed to parse {f}: {e}")

    combined = pd.concat(dfs, ignore_index=True, sort=False)

    # put the most useful columns first; everything else follows in whatever order pandas found it
    front_cols = ["aliquot_barcode", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                  "SVTYPE", "SCTG", "AA_SCTG_TIER", "AA_AMPLICON", "AA_CYCLE", "AA_CYCLECLASS",
                  "GT", "AD", "DP"]
    front_cols = [c for c in front_cols if c in combined.columns]
    other_cols = [c for c in combined.columns if c not in front_cols]
    combined = combined[front_cols + other_cols]

    os.makedirs(os.path.dirname(args.out), exist_ok=True) if os.path.dirname(args.out) else None
    combined.to_csv(args.out, index=False)
    print(f"Wrote combined CSV: {args.out}  ({combined.shape[0]} rows, {combined.shape[1]} columns)")
    print(f"Samples included: {combined['aliquot_barcode'].nunique()}")


if __name__ == "__main__":
    main()

