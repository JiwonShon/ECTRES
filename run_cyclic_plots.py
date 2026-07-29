#!/usr/bin/env python3
"""
run_cyclic_plots.py
────────────────────────────────────────────────────────────────────────────────
aaSuite summary CSV에서 두 가지 조건으로 amplicon_aa_style_v5.R 를 실행한다.

  [Job 1] amplicon_type == 'ecDNA'         → --only_cyclic TRUE  → OUT_DIR/
  [Job 2] OncogenesAmplified 에 ERBB2 포함 → --only_cyclic FALSE → OUT_DIR/ERBB2/

Usage
─────
python3 run_cyclic_plots.py \\
    -i aaSuite_gemline_ms_all_20260430.csv \\
    -s /path/to/amplicon_aa_style_v5.R \\
    -b /mnt/NAS3/.../calls \\
    [-o OUT_DIR] [-j JOBS] [-n TOP_N] [--dry-run]
"""

import argparse
import csv
import logging
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

# ── defaults ──────────────────────────────────────────────────────────────────
DEFAULT_OUT_DIR = (
    "/mnt/NAS3/home/jiwon/ECTRES/summary/aaSuite_germline_ms"
    "/cycle_plot/amplicon_aa_style_v5"
)
DEFAULT_JOBS = 4
DEFAULT_TOP_N = "Inf"

# ── logging ───────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s %(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ── arg parse ─────────────────────────────────────────────────────────────────
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--input",    required=True, metavar="CSV",
                   help="aaSuite summary CSV")
    p.add_argument("-s", "--rscript",  required=True, metavar="R",
                   help="Path to amplicon_aa_style_v5.R")
    p.add_argument("-b", "--base-dir", required=True, metavar="DIR",
                   help="Base directory of aaSuite call results")
    p.add_argument("-o", "--out-dir",  default=DEFAULT_OUT_DIR, metavar="DIR",
                   help=f"Root output directory [default: {DEFAULT_OUT_DIR}]")
    p.add_argument("-j", "--jobs",     default=DEFAULT_JOBS, type=int,
                   help=f"Parallel jobs [default: {DEFAULT_JOBS}]")
    p.add_argument("-n", "--top-n",    default=DEFAULT_TOP_N, metavar="N",
                   help=f"Max cycles (Inf = all) [default: {DEFAULT_TOP_N}]")
    p.add_argument("--dry-run", action="store_true",
                   help="Print commands without executing")
    return p.parse_args()


# ── CSV parsing ───────────────────────────────────────────────────────────────
def load_records(csv_path: str) -> tuple[list[dict], list[dict]]:
    """
    Returns:
        ecdna_records  : amplicon_type == 'ecDNA'
        erbb2_records  : OncogenesAmplified contains 'ERBB2'
    """
    ecdna, erbb2 = [], []
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            barcode = row["amplicon_barcode"].strip()
            # ECTRES-...-amplicon5  →  ECTRES-..._amplicon5
            prefix  = re.sub(r"-amplicon(\d+)$", r"_amplicon\1", barcode)
            rec     = {"prefix": prefix, **row}

            if row.get("amplicon_type", "").strip() == "ecDNA":
                ecdna.append(rec)

            oncogenes = row.get("OncogenesAmplified", "")
            if "ERBB2" in [g.strip() for g in oncogenes.split(",")]:
                erbb2.append(rec)

    return ecdna, erbb2


# ── run one amplicon ──────────────────────────────────────────────────────────
def build_cmd(
    record: dict,
    rscript: str,
    base_dir: str,
    out_dir: Path,
    top_n: str,
    only_cyclic: bool,
) -> list[str]:
    cmd = [
        "Rscript", rscript,
        f"--prefix={record['prefix']}",
        f"--base_dir={base_dir}",
        f"--out_dir={out_dir}",
        f"--top_n={top_n}",
    ]
    if only_cyclic:
        cmd.append("--only_cyclic")   # store_true flag: 값 없이 플래그만
    return cmd


def run_one(record: dict, rscript: str, base_dir: str, out_dir: Path,
            top_n: str, only_cyclic: bool, dry_run: bool) -> dict:
    prefix   = record["prefix"]
    cmd      = build_cmd(record, rscript, base_dir, out_dir, top_n, only_cyclic)
    log_path = out_dir / f"{prefix}.log"

    if dry_run:
        log.info("[DRY-RUN] %s", " ".join(cmd))
        return {"prefix": prefix, "status": "dry-run", "log": None}

    log.info("[START] %s", prefix)
    try:
        with open(log_path, "w") as lf:
            result = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, text=True)
        if result.returncode == 0:
            log.info("[OK   ] %s", prefix)
            return {"prefix": prefix, "status": "ok",    "log": str(log_path)}
        else:
            log.warning("[ERR  ] %s  → see %s", prefix, log_path)
            return {"prefix": prefix, "status": "error", "log": str(log_path)}
    except Exception as exc:
        log.error("[FATAL] %s: %s", prefix, exc)
        return {"prefix": prefix, "status": "fatal", "log": str(log_path)}


# ── batch runner ──────────────────────────────────────────────────────────────
def run_batch(
    label: str,
    records: list[dict],
    rscript: str,
    base_dir: str,
    out_dir: Path,
    top_n: str,
    only_cyclic: bool,
    jobs: int,
    dry_run: bool,
) -> list[dict]:
    out_dir.mkdir(parents=True, exist_ok=True)
    cyclic_str = "TRUE" if only_cyclic else "FALSE"
    log.info("━" * 60)
    log.info("[%s]  %d amplicons  |  only_cyclic=%s  |  out_dir=%s",
             label, len(records), cyclic_str, out_dir)
    log.info("━" * 60)

    results = []
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {
            pool.submit(
                run_one, rec, rscript, base_dir, out_dir, top_n, only_cyclic, dry_run
            ): rec["prefix"]
            for rec in records
        }
        for future in as_completed(futures):
            results.append(future.result())

    ok    = sum(1 for r in results if r["status"] == "ok")
    err   = sum(1 for r in results if r["status"] in ("error", "fatal"))
    log.info("[%s] done  OK=%d  ERR=%d", label, ok, err)
    return results


# ── summary TSV ───────────────────────────────────────────────────────────────
def write_summary(results: list[dict], out_dir: Path, label: str) -> None:
    ts   = datetime.now().strftime("%Y%m%d_%H%M%S")
    path = out_dir / f"run_summary_{label}_{ts}.tsv"
    with open(path, "w") as f:
        f.write("prefix\tstatus\tlog\n")
        for r in results:
            f.write(f"{r['prefix']}\t{r['status']}\t{r.get('log','')}\n")
    log.info("Summary: %s", path)


# ── main ──────────────────────────────────────────────────────────────────────
def main() -> None:
    args = parse_args()

    if not args.dry_run:
        for label, path, check in [
            ("input CSV",  args.input,    Path(args.input).is_file),
            ("Rscript",    args.rscript,  Path(args.rscript).is_file),
            ("base_dir",   args.base_dir, Path(args.base_dir).is_dir),
        ]:
            if not check():
                log.error("%s not found: %s", label, path)
                sys.exit(1)

    root     = Path(args.out_dir)
    erbb2_dir = root / "ERBB2"

    # ── load ──
    ecdna_records, erbb2_records = load_records(args.input)
    log.info("ecDNA amplicons : %d", len(ecdna_records))
    log.info("ERBB2 amplicons : %d", len(erbb2_records))

    if not ecdna_records and not erbb2_records:
        log.warning("Nothing to run.")
        sys.exit(0)

    all_results = {}

    # ── Job 1: ecDNA → only_cyclic TRUE → root out_dir ──
    if ecdna_records:
        res = run_batch(
            label="ecDNA",
            records=ecdna_records,
            rscript=args.rscript,
            base_dir=args.base_dir,
            out_dir=root,
            top_n=args.top_n,
            only_cyclic=True,
            jobs=args.jobs,
            dry_run=args.dry_run,
        )
        all_results["ecDNA"] = res
        if not args.dry_run:
            write_summary(res, root, "ecDNA")

    # ── Job 2: ERBB2 → only_cyclic FALSE → root/ERBB2/ ──
    if erbb2_records:
        res = run_batch(
            label="ERBB2",
            records=erbb2_records,
            rscript=args.rscript,
            base_dir=args.base_dir,
            out_dir=erbb2_dir,
            top_n=args.top_n,
            only_cyclic=False,
            jobs=args.jobs,
            dry_run=args.dry_run,
        )
        all_results["ERBB2"] = res
        if not args.dry_run:
            write_summary(res, erbb2_dir, "ERBB2")

    # ── final summary ──
    print()
    log.info("═" * 60)
    if args.dry_run:
        log.info("  *** DRY-RUN MODE — nothing was actually executed ***")
        log.info("  Remove --dry-run to run for real.")
        log.info("─" * 60)
    for label, res in all_results.items():
        ok  = sum(1 for r in res if r["status"] == "ok")
        err = sum(1 for r in res if r["status"] in ("error", "fatal"))
        dry = sum(1 for r in res if r["status"] == "dry-run")
        if args.dry_run:
            log.info("  %-8s  total=%-4d  DRY-RUN=%d", label, len(res), dry)
        else:
            log.info("  %-8s  total=%-4d  OK=%-4d  ERR=%d", label, len(res), ok, err)
    log.info("═" * 60)

    failed = sum(
        1 for res in all_results.values()
        for r in res if r["status"] in ("error", "fatal")
    )
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()