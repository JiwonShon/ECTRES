import pandas as pd
import glob
import os
import argparse
from datetime import datetime

def parse_params_files(base_dir):
    """
    Recursively find and parse *.params.txt files.
    Extract only the upper summary block starting from 'Gender:'
    and stop before the likelihood table starting with 'n_0' or blank line.
    """
    param_files = glob.glob(os.path.join(base_dir, '**', '*.params.txt'), recursive=True)
    all_results = []

    for file in param_files:
        sample_info = {}

        with open(file, 'r') as f:
            lines = f.readlines()

            for i, line in enumerate(lines):
                if "Gender:" in line:
                    # Sample name is assumed to be the previous line
                    if i > 0:
                        sample_info['Sample'] = lines[i - 1].strip()
                    else:
                        sample_info['Sample'] = os.path.basename(file).replace('.params.txt', '')

                    # Parse key:value lines from Gender: onward
                    for j in range(i, len(lines)):
                        current_line = lines[j].strip()

                        # Stop at blank line or lower likelihood table
                        if not current_line or current_line.startswith('n_0'):
                            break

                        if ":" in current_line:
                            key, val = current_line.split(":", 1)
                            sample_info[key.strip()] = val.strip()
                    break

        if sample_info:
            sample_info['aliquot_barcode'] = sample_info['Sample']
            all_results.append(sample_info)

    df_params = pd.DataFrame(all_results)

    if not df_params.empty:
        cols = ['Sample'] + [c for c in df_params.columns if c != 'Sample']
        df_params = df_params[cols]

    return df_params


def main():
    # 1. Input arguments
    parser = argparse.ArgumentParser(description="Merge ichorCNA output files into single CSVs.")
    parser.add_argument('-i', '--base_dir', required=True, help="Input base directory containing ichorCNA results")
    parser.add_argument('-o', '--save_dir', required=True, help="Directory to save the merged CSV files")
    parser.add_argument('-p', '--project_name', required=True, help="Project name to be used as prefix for output files")

    args = parser.parse_args()

    base_dir = args.base_dir
    save_dir = args.save_dir
    project_name = args.project_name

    current_date = datetime.now().strftime("%Y%m%d")

    os.makedirs(save_dir, exist_ok=True)

    print(f"[*] Starting to merge ichorCNA files for project: {project_name}")
    print(f"[*] Base directory: {base_dir}")
    print(f"[*] Save directory: {save_dir}")

    # 2. Search files
    depth_files = glob.glob(os.path.join(base_dir, '**', '*.correctedDepth.txt'), recursive=True)
    seg_files = glob.glob(os.path.join(base_dir, '**', '*.seg.txt'), recursive=True)
    cna_seg_files = glob.glob(os.path.join(base_dir, '**', '*.cna.seg'), recursive=True)
    param_files = glob.glob(os.path.join(base_dir, '**', '*.params.txt'), recursive=True)

    print(
        f"[-] Found {len(depth_files)} correctedDepth files, "
        f"{len(seg_files)} seg files, "
        f"{len(cna_seg_files)} cna.seg files, "
        f"{len(param_files)} params files."
    )

    # ==========================================
    # [1] correctedDepth.txt merge
    # ==========================================
    list_depth = []
    for f in depth_files:
        sample_id = os.path.basename(f).replace('.correctedDepth.txt', '')
        df = pd.read_csv(f, sep='\t')
        df.insert(0, 'aliquot_barcode', sample_id)
        list_depth.append(df)

    df_depth_merged = pd.concat(list_depth, ignore_index=True) if list_depth else pd.DataFrame()

    # ==========================================
    # [2] seg.txt merge
    # ==========================================
    list_seg = []
    for f in seg_files:
        df = pd.read_csv(f, sep='\t')
        list_seg.append(df)

    df_seg_merged = pd.concat(list_seg, ignore_index=True) if list_seg else pd.DataFrame()

    # ==========================================
    # [3] cna.seg merge
    # ==========================================
    list_cna_seg = []
    for f in cna_seg_files:
        sample_id = os.path.basename(f).replace('.cna.seg', '')
        df = pd.read_csv(f, sep='\t')

        # Remove 'sample_id.' prefix from column names
        new_columns = {col: col.replace(f"{sample_id}.", "") for col in df.columns if sample_id in col}
        df = df.rename(columns=new_columns)

        df.insert(0, 'sample_id', sample_id)
        list_cna_seg.append(df)

    df_cna_seg_merged = pd.concat(list_cna_seg, ignore_index=True) if list_cna_seg else pd.DataFrame()

    # ==========================================
    # [4] params.txt merge
    # ==========================================
    df_params_merged = parse_params_files(base_dir)

    # ==========================================
    # 3. Save CSVs
    # ==========================================
    depth_out_name = os.path.join(save_dir, f"{project_name}_ichorCNA_correctedDepth_{current_date}.csv")
    seg_out_name = os.path.join(save_dir, f"{project_name}_ichorCNA_seg_{current_date}.csv")
    cna_seg_out_name = os.path.join(save_dir, f"{project_name}_ichorCNA_cna_seg_{current_date}.csv")
    params_out_name = os.path.join(save_dir, f"{project_name}_ichorCNA_params_{current_date}.csv")

    if not df_depth_merged.empty:
        df_depth_merged.to_csv(depth_out_name, index=False)
        print(f"[+] Saved: {depth_out_name}")

    if not df_seg_merged.empty:
        df_seg_merged.to_csv(seg_out_name, index=False)
        print(f"[+] Saved: {seg_out_name}")

    if not df_cna_seg_merged.empty:
        df_cna_seg_merged.to_csv(cna_seg_out_name, index=False)
        print(f"[+] Saved: {cna_seg_out_name}")

    if not df_params_merged.empty:
        df_params_merged.to_csv(params_out_name, index=False)
        print(f"[+] Saved: {params_out_name}")

    print("[*] Merge process completed successfully!")


if __name__ == '__main__':
    main()


# Example
# python merge_ichorCNA.py \
#   --base_dir /mnt/NAS3/home/jiwon/HL-NF/scratch/ECTRES/results/ichorCNA/v2/ichorCNA/outlier/ \
#   --save_dir /mnt/NAS3/home/jiwon/ECTRES/summary/ichorCNA/outlier/ \
#   --project_name ECTRES