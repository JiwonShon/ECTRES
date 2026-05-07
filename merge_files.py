import os
import sys
import pandas as pd
from glob import glob


def extract_aliquot_barcode(file_path):
    """
    Extract the aliquot barcode from the sample directory name in the file path.
    Walks up the directory tree to find the sample-level folder
    (e.g., ECTRES-ECGI1-0001-TPX-A01-WGS-1ST985).
    
    Identifies sample directories by the pattern: 6+ dash-separated parts
    where the first part looks like a project code (all uppercase letters).
    """
    parts = file_path.split(os.sep)
    for part in reversed(parts):
        segments = part.split("-")
        if len(segments) >= 6 and segments[0].isupper() and segments[0].isalpha():
            return part
    return "Unknown"


def merge_files(input_folder, file_pattern, file_format, output_dir):
    """
    Parameters:
        input_folder (str): Root folder containing sample subdirectories
        file_pattern (str): File suffix pattern to match (e.g., 'purity')
        file_format (str): File extension without dot (e.g., 'tsv', 'csv', 'txt')
        output_dir (str): Directory to save the merged output and file list
    """
    print("input_folder:", input_folder)
    print("file_pattern:", file_pattern)
    print("file_format:", file_format)
    print("output_dir:", output_dir)

    file_list = glob(
        os.path.join(input_folder, "**", f"*{file_pattern}.{file_format}"),
        recursive=True
    )

    if not file_list:
        print("No files found matching the given pattern.")
        return

    print(f"Found {len(file_list)} file(s):")
    for f in file_list:
        print(" ", f)

    merged_data = []

    for file_path in file_list:
        try:
            if file_format == "tsv":
                df = pd.read_csv(file_path, sep='\t')
            elif file_format == "csv":
                df = pd.read_csv(file_path)
            elif file_format == "txt":
                with open(file_path, "r", encoding="utf-8") as f:
                    lines = f.readlines()
                df = pd.DataFrame({"content": lines})
            else:
                print(f"Unsupported file format: {file_format}")
                continue

            aliquot_barcode = extract_aliquot_barcode(file_path)
            df.insert(0, "aliquot_barcode", aliquot_barcode)

            merged_data.append(df)

        except Exception as e:
            print(f"Error reading file {file_path}: {e}")

    if not merged_data:
        print("No data to merge.")
        return

    merged_df = pd.concat(merged_data, ignore_index=True)

    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, f"combined.{file_pattern}.csv")
    merged_df.to_csv(output_file, index=False)
    print(f"\nMerged file saved: {output_file}")

    file_list_path = os.path.join(output_dir, "merged_file_list.txt")
    with open(file_list_path, "w") as f:
        for file in file_list:
            f.write(file + "\n")
    print(f"File list saved: {file_list_path}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python merge_files.py <input_folder> <file_pattern> <file_format> <output_dir>")
        print("Example: python merge_files.py /path/to/samples purity tsv /path/to/output")
        sys.exit(1)

    input_folder = sys.argv[1]
    file_pattern = sys.argv[2]
    file_format = sys.argv[3]
    output_dir = sys.argv[4]

    merge_files(input_folder, file_pattern, file_format, output_dir)
