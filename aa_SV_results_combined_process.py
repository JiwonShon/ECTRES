import pandas as pd
import glob
import os
import argparse
from datetime import datetime

def main():
    # 1. 명령어 인수 설정
    parser = argparse.ArgumentParser(description="Merge AA SV_summary files into a single CSV.")
    parser.add_argument('-i', '--base_dir', required=True, help="Input base directory containing AA results")
    parser.add_argument('-o', '--save_dir', required=True, help="Directory to save the merged CSV files")
    parser.add_argument('-p', '--project_name', required=True, help="Project name to be used as prefix for output files")
    
    args = parser.parse_args()
    
    base_dir = args.base_dir
    save_dir = args.save_dir
    project_name = args.project_name
    
    # 오늘 날짜 추출
    current_date = datetime.now().strftime("%Y%m%d")
    
    os.makedirs(save_dir, exist_ok=True)

    print(f"[*] Starting to merge SV_summary files for project: {project_name}")
    print(f"[*] Base directory: {base_dir}")
    print(f"[*] Save directory: {save_dir}")

    # 2. SV_summary.tsv 파일 검색
    sv_files = glob.glob(os.path.join(base_dir, '**', '*_SV_summary.tsv'), recursive=True)
    print(f"[-] Found {len(sv_files)} SV_summary.tsv files.")

    # ==========================================
    # [1] 파일 파싱 및 병합 (Long Format)
    # ==========================================
    list_sv = []
    for f in sv_files:
        # 파일명 추출 (예: ECTRES-EFM19-0001-TPX-A01-WGS-HCS4K_amplicon6_SV_summary.tsv)
        base_name = os.path.basename(f)
        
        # 1. amplicon_barcode 추출: '_SV_summary.tsv' 제거
        # 결과: ECTRES-EFM19-0001-TPX-A01-WGS-HCS4K_amplicon6
        amplicon_barcode = base_name.replace('_SV_summary.tsv', '')
        
        # 2. aa_barcode와 amplicon_number 분리
        # 오른쪽(뒤)에서부터 첫 번째 '_'를 기준으로 쪼갬 (rsplit)
        if '_' in amplicon_barcode:
            aa_barcode, amplicon_number = amplicon_barcode.rsplit('_', 1)
        else:
            aa_barcode = amplicon_barcode
            amplicon_number = 'unknown'

        # TSV 파일 읽기
        df = pd.read_csv(f, sep='\t')
        
        # 새로운 컬럼들을 맨 앞에 순서대로 삽입 (index 0)
        df.insert(0, 'amplicon_number', amplicon_number)
        df.insert(0, 'aa_barcode', aa_barcode)
        df.insert(0, 'amplicon_barcode', amplicon_barcode)
        
        list_sv.append(df)

    # 전체 리스트 병합
    if list_sv:
        df_sv_merged = pd.concat(list_sv, ignore_index=True)
    else:
        df_sv_merged = pd.DataFrame()

    # ==========================================
    # [2] CSV로 저장
    # ==========================================
    sv_out_name = os.path.join(save_dir, f"{project_name}-aaSuite_germline_ms-v1.3.8-GRCh37-minCN4.5-cnsizeMin50000-10X_SV_summary_{current_date}.csv")

    if not df_sv_merged.empty:
        df_sv_merged.to_csv(sv_out_name, index=False)
        print(f"[+] Saved: {sv_out_name}")
    else:
        print("[!] No data to save. Please check the input directory.")

    print("[*] Merge process completed successfully!")

if __name__ == '__main__':
    main()

# python aa_SV_results_combined_process.py \
#   --base_dir /mnt/NAS3/home/jiwon/HL-NF/scratch/ECTRES/results/aaSuite_germline_ms/v1.3.8/GRCh37/minCN4.5/cnsizeMin50000/10X/calls/ \
#   --save_dir /mnt/NAS3/home/jiwon/ECTRES/summary/aaSuite_germline_ms/10X/ \
#   --project_name ECTRES


