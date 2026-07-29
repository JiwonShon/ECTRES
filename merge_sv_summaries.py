#!/usr/bin/env python3
import os
import re
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Merge AA SV summary TSV files")
parser.add_argument("--base_dir", required=True, help="AA calls 최상위 디렉토리 경로")
parser.add_argument("--output_dir", required=True, help="결과 파일 저장 디렉토리 경로")
args = parser.parse_args()

BASE_DIR = args.base_dir
OUTPUT_DIR = args.output_dir
os.makedirs(OUTPUT_DIR, exist_ok=True)

# SV summary 파일 패턴
SV_PATTERN = re.compile(r"^(.+)_amplicon(\d+)_SV_summary\.tsv$")

dfs = []
skipped = []

for barcode in os.listdir(BASE_DIR):
    barcode_dir = os.path.join(BASE_DIR, barcode)
    if not os.path.isdir(barcode_dir):
        continue

    # 하위 디렉토리 재귀 탐색
    for root, dirs, files in os.walk(barcode_dir):
        for fname in files:
            m = SV_PATTERN.match(fname)
            if not m:
                continue

            aa_barcode = m.group(1)
            amplicon_number = f"amplicon{m.group(2)}"
            fpath = os.path.join(root, fname)

            try:
                df = pd.read_csv(fpath, sep="\t")
            except Exception as e:
                print(f"[WARN] 읽기 실패: {fpath} ({e})")
                continue

            # 데이터 없이 컬럼만 있는 경우 skip
            if df.empty:
                skipped.append(fpath)
                print(f"[SKIP] 데이터 없음: {fpath}")
                continue

            df.insert(0, "aa_barcode", aa_barcode)
            df.insert(1, "amplicon_number", amplicon_number)
            dfs.append(df)

if not dfs:
    print("합칠 파일이 없어요.")
else:
    merged = pd.concat(dfs, ignore_index=True)
    out_path = os.path.join(OUTPUT_DIR, "merged_SV_summary_all.tsv")
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"\n✅ 완료: {len(dfs)}개 파일 합침")
    print(f"✅ 총 {len(merged)}개 SV rows")
    print(f"✅ 저장 위치: {out_path}")
    if skipped:
        print(f"\n⚠️  빈 파일 {len(skipped)}개 skip됨:")
        for s in skipped:
            print(f"   {s}")


#python3 merge_sv_summaries.py \
#  --base_dir /mnt/NAS3/home/jiwon/HL-NF/scratch/ECTRES/results/aaSuite_germline_ms/v1.3.8/GRCh37/minCN4.5/cnsizeMin50000/10X/calls \
#  --output_dir /mnt/NAS3/home/jiwon/HL-NF/scratch/ECTRES/results/merged_output
