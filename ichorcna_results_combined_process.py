import pandas as pd
import glob
import os
import argparse
from datetime import datetime

def main():
    # 1. 명령어 인수 설정 (Input arguments)
    parser = argparse.ArgumentParser(description="Merge ichorCNA output files into single CSVs.")
    parser.add_argument('-i', '--base_dir', required=True, help="Input base directory containing ichorCNA results")
    parser.add_argument('-o', '--save_dir', required=True, help="Directory to save the merged CSV files")
    parser.add_argument('-p', '--project_name', required=True, help="Project name to be used as prefix for output files")
    
    args = parser.parse_args()
    
    base_dir = args.base_dir
    save_dir = args.save_dir
    project_name = args.project_name
    
    # 오늘 날짜 추출 (YYYYMMDD 형식)
    current_date = datetime.now().strftime("%Y%m%d")
    
    # 저장할 디렉토리가 없으면 생성
    os.makedirs(save_dir, exist_ok=True)

    print(f"[*] Starting to merge ichorCNA files for project: {project_name}")
    print(f"[*] Base directory: {base_dir}")
    print(f"[*] Save directory: {save_dir}")

    # 2. 파일 검색
    depth_files = glob.glob(os.path.join(base_dir, '**', '*.correctedDepth.txt'), recursive=True)
    seg_files = glob.glob(os.path.join(base_dir, '**', '*.seg.txt'), recursive=True)
    cna_seg_files = glob.glob(os.path.join(base_dir, '**', '*.cna.seg'), recursive=True)

    print(f"[-] Found {len(depth_files)} correctedDepth files, {len(seg_files)} seg files, {len(cna_seg_files)} cna.seg files.")

    # ==========================================
    # [1] correctedDepth.txt 병합 (Long Format으로 변경)
    # ==========================================
    list_depth = []
    for f in depth_files:
        sample_id = os.path.basename(f).replace('.correctedDepth.txt', '')
        df = pd.read_csv(f, sep='\t')
        df.insert(0, 'aliquot_barcode', sample_id)
        list_depth.append(df)

    df_depth_merged = pd.concat(list_depth, ignore_index=True) if list_depth else pd.DataFrame()
    
    # ==========================================
    # [2] seg.txt 병합 (Long Format)
    # ==========================================
    list_seg = []
    for f in seg_files:
        df = pd.read_csv(f, sep='\t')
        list_seg.append(df)
    
    df_seg_merged = pd.concat(list_seg, ignore_index=True) if list_seg else pd.DataFrame()

    # ==========================================
    # [3] cna.seg 병합 (Long Format, 컬럼명 정리)
    # ==========================================
    list_cna_seg = []
    for f in cna_seg_files:
        sample_id = os.path.basename(f).replace('.cna.seg', '')
        df = pd.read_csv(f, sep='\t')
        
        # 컬럼명에서 'sample_id.' 부분 제거
        new_columns = {col: col.replace(f"{sample_id}.", "") for col in df.columns if sample_id in col}
        df = df.rename(columns=new_columns)
        
        # 어느 샘플인지 식별하기 위해 맨 앞에 sample_id 컬럼 추가
        df.insert(0, 'sample_id', sample_id)
        list_cna_seg.append(df)

    df_cna_seg_merged = pd.concat(list_cna_seg, ignore_index=True) if list_cna_seg else pd.DataFrame()

    # ==========================================
    # 3. CSV로 저장 (구분자 콤마, 프로젝트명 & 날짜 포함)
    # ==========================================
    # 파일명 생성
    depth_out_name = os.path.join(save_dir, f"{project_name}_ichorCNA_correctedDepth_{current_date}.csv")
    seg_out_name = os.path.join(save_dir, f"{project_name}_ichorCNA_seg_{current_date}.csv")
    cna_seg_out_name = os.path.join(save_dir, f"{project_name}_ichorCNA_cna_seg_{current_date}.csv")

    # CSV 저장 (sep='\t' 제거하여 기본값인 콤마(,) 사용)
    if not df_depth_merged.empty:
        df_depth_merged.to_csv(depth_out_name, index=False)
        print(f"[+] Saved: {depth_out_name}")
        
    if not df_seg_merged.empty:
        df_seg_merged.to_csv(seg_out_name, index=False)
        print(f"[+] Saved: {seg_out_name}")
        
    if not df_cna_seg_merged.empty:
        df_cna_seg_merged.to_csv(cna_seg_out_name, index=False)
        print(f"[+] Saved: {cna_seg_out_name}")

    print("[*] Merge process completed successfully!")

if __name__ == '__main__':
    main()


# python merge_ichorCNA.py \
#   --base_dir /mnt/NAS3/home/jiwon/HL-NF/scratch/ECTRES/results/ichorCNA/v2/ichorCNA/ \
#   --save_dir /mnt/NAS3/home/jiwon/ECTRES/summary/ichorCNA/ \
#   --project_name ECTRES
