#!/usr/bin/env python3
"""
auto_cycleviz.py

For a given aliquot_barcode, this scans AmpliconClassifier's
"_classification_bed_files" directory for every "*_ecDNA_N_intervals.bed"
file (across all amplicons, and all ecDNA within an amplicon), then
figures out WHICH "Cycle=" entries in the matching AA cycles.txt actually
correspond to that ecDNA, by comparing genomic overlap between each
cycle's segments and the bed file's intervals. AA/AC often reports
several alternative cyclic decompositions for the same ecDNA locus
(e.g. Cycle=1, Cycle=2, Cycle=3, all tagged ecDNA-like with similar
Copy_count but different segment combinations) - the ecDNA_N bed file
is the union of all of them. This script plots EVERY matching
decomposition (not just one "best" pick) with the correct --cycle
number filled in automatically, so you can compare them. Use
--best-only if you only want the single highest-overlap decomposition.

Why this is needed: AmpliconClassifier's bed filenames only tell you
"amplicon25_ecDNA_1", they don't tell you which numbered Cycle(s) in
amplicon25_cycles.txt correspond to it. This script figures that out
by matching genomic coordinates instead of you reading the cycles file
by hand each time.

Usage:
  python3 auto_cycleviz.py -b <barcode> [--dry-run]

  # sanity-check the matches first, without running CycleViz:
  python3 auto_cycleviz.py -b ECTRES-ECGI1-0001-TPX-A07-WGS-7ZJ410 --dry-run

  # restrict to one amplicon:
  python3 auto_cycleviz.py -b ECTRES-H2170-0001-TPX-A01-WGS-3YV111 --amplicon 12

Assumes the standard AmpliconSuite/AmpliconClassifier directory layout:
  <data_root>/<barcode>_output/<barcode>_classification/<barcode>_classification_bed_files/
  <data_root>/<barcode>_output/<barcode>_AA_results/
"""

import argparse
import glob
import os
import re
import subprocess
import sys
from collections import defaultdict


def parse_cycles_file(path):
    """Parse an AA/AC cycles.txt file.
    Returns (segments: dict[int -> (chrom, start, end)], cycles: list[dict])
    """
    segments = {}
    cycles = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("Segment"):
                parts = line.split()
                # Segment <id> <chrom> <start> <end>
                seg_id = int(parts[1])
                chrom = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                segments[seg_id] = (chrom, start, end)
            elif line.startswith("Cycle="):
                fields = dict(kv.split("=", 1) for kv in line.split(";") if "=" in kv)
                segments_raw = fields.get("Segments", "")
                if "IsCyclicPath" in fields:
                    is_cyclic = fields["IsCyclicPath"] == "True"
                else:
                    # Plain AA cycles.txt has no IsCyclicPath/CycleClass tags (those are
                    # only added by AmpliconClassifier's annotated_cycles_files). Fall back
                    # to AA's own convention: linear paths start/end at the reference
                    # segment "0"; circular/cyclic paths never touch segment 0.
                    toks = [t for t in segments_raw.split(",") if t]
                    is_cyclic = toks and not any(int(t[:-1]) == 0 for t in toks)
                cyc = {
                    "id": int(fields["Cycle"]),
                    "copy_count": float(fields.get("Copy_count", 0)),
                    "length": int(fields.get("Length", 0)),
                    "is_cyclic": bool(is_cyclic),
                    "class": fields.get("CycleClass", ""),
                    "segments_raw": segments_raw,
                }
                cycles.append(cyc)
    return segments, cycles


def cycle_intervals(cycle, segments):
    """Genomic intervals visited by a cycle (skips the 0/telomere marker)."""
    ivals = []
    for tok in cycle["segments_raw"].split(","):
        tok = tok.strip()
        if not tok:
            continue
        seg_id = int(tok[:-1])
        if seg_id == 0:
            continue
        if seg_id in segments:
            ivals.append(segments[seg_id])
    return ivals


def read_bed(path):
    ivals = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            f = line.split()
            ivals.append((f[0], int(f[1]), int(f[2])))
    return ivals


def overlap_bp(a, b):
    if a[0] != b[0]:
        return 0
    return max(0, min(a[2], b[2]) - max(a[1], b[1]))


def merge_intervals(ivals):
    """Merge overlapping/touching (chrom, start, end) intervals within each chrom."""
    by_chrom = defaultdict(list)
    for c, s, e in ivals:
        by_chrom[c].append((s, e))

    merged = []
    for c, lst in by_chrom.items():
        lst.sort()
        cur_s, cur_e = lst[0]
        for s, e in lst[1:]:
            if s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                merged.append((c, cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((c, cur_s, cur_e))
    return merged


def jaccard_score(cyc_ivals, bed_ivals, padding=10000):
    """
    Jaccard-style overlap between a cycle's genomic footprint and a bed
    file's footprint: intersection / union. This rewards cycles that
    cover MOST of the bed region (recall), while still penalizing
    cycles that mostly lie outside of it (precision) - a plain
    "is the cycle a subset of bed" check isn't enough, since an
    ecDNA_N bed file is often the UNION of segments from several
    alternative AA cycle decompositions of the same ecDNA, so a small
    single-segment cycle can be 100% "contained" in the bed while only
    explaining a small fraction of the actual structure.
    """
    padded_cyc = merge_intervals([(c, max(0, s - padding), e + padding) for c, s, e in cyc_ivals])
    padded_bed = merge_intervals([(c, max(0, s - padding), e + padding) for c, s, e in bed_ivals])

    cyc_bp = sum(e - s for _, s, e in padded_cyc)
    bed_bp = sum(e - s for _, s, e in padded_bed)

    inter_bp = 0
    for civ in padded_cyc:
        for biv in padded_bed:
            inter_bp += overlap_bp(civ, biv)

    union_bp = cyc_bp + bed_bp - inter_bp
    if union_bp == 0:
        return 0.0
    return inter_bp / union_bp


def find_matching_cycles(bed_ivals, cycles, segments, min_overlap=0.05):
    """
    Return every cyclic candidate that overlaps the ecDNA bed region above
    min_overlap, as a list of (jaccard_score, cycle_dict), best match first.

    AA/AC frequently reports several alternative cyclic decompositions for
    the same ecDNA locus (near-identical Copy_count, different segment
    combinations e.g. Cycle=1, Cycle=2, Cycle=3) - the ecDNA_N bed file is
    the UNION of segments across all of them. Rather than silently picking
    one "winner", we surface every plausible decomposition so all of them
    can be plotted and compared.
    """
    candidates = []
    for cyc in cycles:
        if not cyc["is_cyclic"]:
            continue
        civals = cycle_intervals(cyc, segments)
        if not civals:
            continue
        score = jaccard_score(civals, bed_ivals)
        if score >= min_overlap:
            candidates.append((score, cyc))

    candidates.sort(key=lambda x: -x[0])
    return candidates


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-b", "--barcode", required=True)
    ap.add_argument("-n", "--out-name", default=None,
                     help="label to use in OUTPUT filenames instead of the barcode (e.g. sample_id). "
                          "Input paths are still located using --barcode.")
    ap.add_argument("-d", "--data-root", default=os.environ.get("CALLS_ROOT", "/mnt/AA_calls"),
                     help="dir containing <barcode>/<barcode>_output/ (default: $CALLS_ROOT or /mnt/AA_calls)")
    ap.add_argument("-o", "--results-root", default="/mnt/ECTRES/results/CycleViz")
    ap.add_argument("-r", "--ref", default="GRCh37")
    ap.add_argument("-f", "--gene-fontsize", default="10")
    ap.add_argument("-x", "--extra-args", default="--noPDF", help="extra args passed through to CycleViz.py")
    ap.add_argument("--min-overlap", type=float, default=0.05,
                     help="minimum Jaccard overlap (intersection/union) between a cycle's genomic "
                          "footprint and the ecDNA bed file's footprint to be included as a match "
                          "(default 0.05 - low, since we want every plausible decomposition, not "
                          "just the single best one)")
    ap.add_argument("--best-only", action="store_true",
                     help="only run CycleViz for the single highest-scoring cycle per ecDNA feature, "
                          "instead of every cyclic decomposition that overlaps it")
    ap.add_argument("--cv-src", default=os.environ.get("CV_SRC"), help="path to CycleViz repo (default: $CV_SRC)")
    ap.add_argument("--dry-run", action="store_true", help="only print the amplicon/ecDNA -> cycle matches, don't run CycleViz")
    ap.add_argument("--amplicon", help="restrict to a single amplicon number, e.g. 25")
    args = ap.parse_args()

    if not args.cv_src and not args.dry_run:
        sys.exit("CV_SRC is not set. `export CV_SRC=/path/to/CycleViz` or pass --cv-src, or use --dry-run.")

    out_label = args.out_name if args.out_name else args.barcode

    out_base = args.barcode + "_output"

    # Support both layouts:
    #   <data_root>/<barcode>/<barcode>_output/...   (e.g. the AA "calls" dir mounted as-is)
    #   <data_root>/<barcode>_output/...              (e.g. a single sample copied into data/)
    candidate_roots = [
        os.path.join(args.data_root, args.barcode),
        args.data_root,
    ]
    class_dir = aa_dir = None
    for root in candidate_roots:
        cd = os.path.join(root, out_base, args.barcode + "_classification",
                           args.barcode + "_classification_bed_files")
        if os.path.isdir(cd):
            class_dir = cd
            aa_dir = os.path.join(root, out_base, args.barcode + "_AA_results")
            break

    if class_dir is None:
        tried = "\n  ".join(os.path.join(r, out_base, args.barcode + "_classification",
                                          args.barcode + "_classification_bed_files") for r in candidate_roots)
        sys.exit(f"classification bed dir not found. Tried:\n  {tried}")

    pattern = os.path.join(class_dir, f"{args.barcode}_amplicon*_ecDNA_*_intervals.bed")
    bed_files = sorted(glob.glob(pattern))
    if not bed_files:
        sys.exit(f"no ecDNA interval bed files found under {class_dir}")

    rx = re.compile(re.escape(args.barcode) + r"_amplicon(\d+)_ecDNA_(\d+)_intervals\.bed$")

    cycles_cache = {}
    if not args.dry_run:
        os.makedirs(args.results_root, exist_ok=True)

    print(f"{'amplicon':>8} {'ecDNA':>5} {'cycle':>5} {'copy_count':>10} {'jaccard':>8}  status")
    for bed_path in bed_files:
        m = rx.search(os.path.basename(bed_path))
        if not m:
            continue
        amp_n, ec_n = m.group(1), m.group(2)
        if args.amplicon and amp_n != args.amplicon:
            continue

        cycles_path = os.path.join(aa_dir, f"{args.barcode}_amplicon{amp_n}_cycles.txt")
        graph_path = os.path.join(aa_dir, f"{args.barcode}_amplicon{amp_n}_graph.txt")

        if not os.path.exists(cycles_path):
            print(f"{amp_n:>8} {ec_n:>5} {'--':>5} {'--':>10} {'--':>8}  cycles.txt missing: {cycles_path}")
            continue

        if cycles_path not in cycles_cache:
            cycles_cache[cycles_path] = parse_cycles_file(cycles_path)
        segments, cycles = cycles_cache[cycles_path]

        bed_ivals = read_bed(bed_path)
        matches = find_matching_cycles(bed_ivals, cycles, segments, args.min_overlap)

        if not matches:
            print(f"{amp_n:>8} {ec_n:>5} {'--':>5} {'--':>10} {'--':>8}  no matching cyclic path found "
                  f"in {cycles_path} (try a lower --min-overlap)")
            continue

        if args.best_only:
            matches = matches[:1]

        for score, cyc in matches:
            print(f"{amp_n:>8} {ec_n:>5} {cyc['id']:>5} {cyc['copy_count']:>10.2f} {score:>8.2f}  ok")

            if args.dry_run:
                continue

            out_prefix = os.path.join(args.results_root,
                                       f"{out_label}_amplicon{amp_n}_ecDNA{ec_n}_cycle{cyc['id']}")
            cmd = [
                os.path.join(args.cv_src, "CycleViz.py"),
                "--structure_bed", bed_path,
                "--graph", graph_path,
                "--cycles", cycles_path,
                "--ref", args.ref,
                "-o", out_prefix,
                "--gene_fontsize", str(args.gene_fontsize),
                "--cycle", str(cyc["id"]),
            ] + args.extra_args.split()

            print("  -> " + " ".join(cmd))
            subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()