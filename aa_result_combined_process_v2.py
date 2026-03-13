"""
File_name: aa_result_combined_process_v2.py

Project:aaSuite_germline_ss result summary 
Author: Jiwon Shon
Created: April 9, 2024
Edit: Sep 12, 2024
"""
# Import libries
import sys
import numpy as np
import pandas as pd
import os
from datetime import datetime
import re
import glob

if len(sys.argv) != 7:
    print("Usage: python3 aa_result_combined_process_v2.py <projectName> <minCN> <cnzise> <downsamples> <NAS> <step>")
    print("for exmple, > python3 aa_result_combined_process_v2.py SMCAD2_v2 minCN4 cnsizeMin5000 10X NAS4 aaSuite_somatic_ms")
    sys.exit(1)

# /mnt/NAS3/home/mary/HL-NF/scratch/ECTRES/results/aaSuite_germline_ms/v1.3.8/GRCh37/minCN4.5/cnsizeMin50000/-1X/calls    
#     python3 aa_result_combined_process_v2.py ECTRES minCN4.5 cnsizeMin50000 -1X NAS3 aaSuite_germline_ms

# /mnt/NAS3/home/jiwon/HL-NF/scratch/ECTRES/results/aaSuite_germline_ms/v1.3.8/GRCh37/minCN4.5/cnsizeMin50000/-1X/calls/
#     python3 aa_result_combined_process_v2.py ECTRES minCN4.5 cnsizeMin50000 -1X NAS3 aaSuite_germline_ms
#     python3 aa_result_combined_process_v2.py ECTRES minCN4.5 cnsizeMin50000 1X NAS3 aaSuite_germline_ms
#     python3 aa_result_combined_process_v2.py ECTRES minCN4.5 cnsizeMin50000 10X NAS3 aaSuite_germline_ms

projectName = sys.argv[1]  
minCN = sys.argv[2]  
cnzise = sys.argv[3]  
downsamples = sys.argv[4]  
NAS = sys.argv[5]  
step = sys.argv[6]

base_path = f'/mnt/{NAS}/home/jiwon/HL-NF/scratch/{projectName}/results/{step}/v1.3.8/GRCh37/{minCN}/{cnzise}/{downsamples}/calls'
save_path = f'/mnt/NAS3/home/jiwon/{projectName}/summary/{step}/{downsamples}/'
os.makedirs(save_path, exist_ok=True)

# ============================================================ #
# ============================================================ #

os.makedirs(save_path, exist_ok=True)
# Setup log file
log_date = datetime.now().strftime("%Y%m%d")
log_filename = f"log_{log_date}.txt"
log_filepath = os.path.join(save_path, log_filename)
# Open log file
log_file = open(log_filepath, 'w')
# Redirect print statements to log file
def log_print(*args, **kwargs):
    print(*args, **kwargs, file=log_file)
log_print("File_name: aa_result_combined_process.py")
log_print("Project: aaSuite_germline_ss result summary")
log_print("Author: Jiwon Shon")
log_print("Created: April 9, 2024")
log_print("============================================================")
log_print()

#
current_directory = os.getcwd()
print("Current Directory:", current_directory)
#
manipulate_date = datetime.now().strftime("%Y%m%d")
print('Date of today:', manipulate_date)
# 
path_parts = base_path.split('/')
relevant_parts = path_parts[-7:-1]
print(relevant_parts) 
relevant_parts.insert(0, path_parts[-9])
print(relevant_parts)
formatted_string = '-'.join(relevant_parts)

# Generate a list of directories in the base path
sample_dirs = [name for name in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, name))]
print(sample_dirs)

# ================= 1. combined_aa_amplicons.csv 
process='combined_aa_amplicons'
all_data = []
for sample_dir in sample_dirs:
    summary_file_path = os.path.join(base_path, sample_dir, f"{sample_dir}_output", f"{sample_dir}_AA_results", f"{sample_dir}_summary.txt")
    
    if os.path.isfile(summary_file_path):
        with open(summary_file_path, 'r') as file:
            for line in file:
                key_value = line.strip().split(' = ')
                if ' ' in key_value[0]:
                    amplicon_part, intervals_part = key_value[0].split(' ', 1)
                    intervals_key = intervals_part.replace("#", "N_")
                else:
                    continue

                if '-----' in line:
                    pass
                elif 'AmpliconID' in line:
                    current_amplicon = {'amplicon_barcode':sample_dir+'-amplicon' + key_value[1],'aa_barcode': sample_dir, intervals_key: key_value[1], 'amplicon_index': 'amplicon' + key_value[1], 'aa_summary_file_path': summary_file_path}
                    all_data.append(current_amplicon)
                else:
                    current_amplicon[intervals_key] = key_value[1]

summary_df = pd.DataFrame(all_data)

final_string = f"{formatted_string}_{manipulate_date}_{process}.csv"
print("save file path: ", save_path+final_string)
summary_df.to_csv(save_path + final_string, index=False)
print(f"{process} file created successfully.")

# ================= 2. combined_amplicon_classification_profiles.csv
process='combined_amplicon_classification_profiles'
all_dfs = []
for sample_dir in sample_dirs:
    print("Doing",sample_dir)
    tsv_file_path = os.path.join(base_path, sample_dir, f"{sample_dir}_output", f"{sample_dir}_classification", f"{sample_dir}_amplicon_classification_profiles.tsv")
    print(tsv_file_path)
    if os.path.isfile(tsv_file_path):
        df = pd.read_csv(tsv_file_path, sep='\t')
        df['amplicon_barcode'] = df['sample_name'] + '-' + df['amplicon_number']
        df['aa_barcode'] = df['sample_name']

        all_dfs.append(df)

final_df = pd.concat(all_dfs, ignore_index=True)
final_df=final_df[['amplicon_barcode', 'aa_barcode', 'sample_name', 'amplicon_number', 'amplicon_decomposition_class',
       'ecDNA+', 'BFB+', 'ecDNA_amplicons']]

final_string_2 = f"{formatted_string}_{manipulate_date}_{process}.csv"
print("save file path: ", save_path+final_string_2)
final_df.to_csv(save_path + final_string_2, index=False)
print(f"{process} file created successfully.")

# ================= 3. combined_ecDNA_counts.csv
process='combined_ecDNA_counts'
all_dfs3 = []
for sample_dir in sample_dirs:
    print("Doing",sample_dir)
    tsv_file_path = os.path.join(base_path, sample_dir, f"{sample_dir}_output", f"{sample_dir}_classification", f"{sample_dir}_ecDNA_counts.tsv")
    if os.path.isfile(tsv_file_path):
        df = pd.read_csv(tsv_file_path, sep='\t')
        if '#sample' in df.columns:
            df = df.rename(columns={'#sample': 'aa_barcode'})

        all_dfs3.append(df)
    else:
        print('File not found:', tsv_file_path)

final_df3 = pd.concat(all_dfs3, ignore_index=True)

final_string_3 = f"{formatted_string}_{manipulate_date}_{process}.csv"
print("save file path: ", save_path+final_string_3)
final_df3.to_csv(save_path + final_string_3, index=False)
print(f"{process} file created successfully.")

# ================= 4. combined_feature_basic_properties.csv
process='combined_feature_basic_properties'
all_dfs4 = []
for sample_dir in sample_dirs:
    print("Doing",sample_dir)
    tsv_file_path = os.path.join(base_path, sample_dir, f"{sample_dir}_output", f"{sample_dir}_classification", f"{sample_dir}_feature_basic_properties.tsv")
    if os.path.isfile(tsv_file_path):
        df = pd.read_csv(tsv_file_path, sep='\t')
        df['aa_barcode'] = df['feature_ID'].str.split('_', n=1).str[0]

        all_dfs4.append(df)
    else:
        print('File not found:', tsv_file_path)

final_df4 = pd.concat(all_dfs4, ignore_index=True)
final_df4=final_df4[[ 'aa_barcode', 'feature_ID', 'captured_region_size_bp', 'median_feature_CN',
       'max_feature_CN', 'borderline_flag']]
final_df4['borderline_flag'].fillna('None', inplace=True)

final_string_4 = f"{formatted_string}_{manipulate_date}_{process}.csv"
print("save file path: ", save_path+final_string_4)
final_df4.to_csv(save_path + final_string_4, index=False)
print(f"{process} file created successfully.")

# ================= 5. combined_feature_entropy.csv
process='combined_feature_entropy'
all_dfs5 = []
for sample_dir in sample_dirs:
    print("Doing",sample_dir)
    tsv_file_path = os.path.join(base_path, sample_dir, f"{sample_dir}_output", f"{sample_dir}_classification", f"{sample_dir}_feature_entropy.tsv")
    if os.path.isfile(tsv_file_path):
        df = pd.read_csv(tsv_file_path, sep='\t')
        df['feature_barcode'] = df['sample'] + '-' + df['amplicon'] + '-' + df['feature']
        df['amplicon_barcode'] = df['sample'] + '-' + df['amplicon']
        df['aa_barcode'] = df['sample']

        all_dfs5.append(df)
    else:
        print('File not found:', tsv_file_path)

final_df5 = pd.concat(all_dfs5, ignore_index=True)
final_df5=final_df5[['feature_barcode', 'amplicon_barcode', 'aa_barcode','sample', 'amplicon', 'feature', 'total_feature_entropy',
       'decomp_entropy', 'Amp_nseg_entropy']]

final_string_5 = f"{formatted_string}_{manipulate_date}_{process}.csv"
print("save file path: ", save_path+final_string_5)
final_df5.to_csv(save_path + final_string_5, index=False)
print(f"{process} file created successfully.")

# ================= 6. combined_gene_list.csv
process='combined_gene_list'
all_dfs6 = []
for sample_dir in sample_dirs:
    print("Doing",sample_dir)
    tsv_file_path = os.path.join(base_path, sample_dir, f"{sample_dir}_output", f"{sample_dir}_classification", f"{sample_dir}_gene_list.tsv")
    if os.path.isfile(tsv_file_path):
        df = pd.read_csv(tsv_file_path, sep='\t')
        df['aa_barcode'] = df['sample_name']

        all_dfs6.append(df)
    else:
        print('File not found:', tsv_file_path)

final_df6 = pd.concat(all_dfs6, ignore_index=True)
final_df6=final_df6[[ 'aa_barcode', 'sample_name', 'amplicon_number', 'feature', 'gene', 'gene_cn',
       'truncated', 'is_canonical_oncogene']]
final_df6['truncated'].fillna('None', inplace=True)

final_string_6 = f"{formatted_string}_{manipulate_date}_{process}.csv"
print("save file path: ", save_path+final_string_6)
final_df6.to_csv(save_path + final_string_6, index=False)
print(f"{process} file created successfully.")

# ================= 7. annotated_cycles_files.csv
process='annotated_cycles_files'
import glob

# for sample_dir in sample_dirs:
#     print("Processing",sample_dir)
#     # Construct the file path for the amplicon classification profiles TSV file
#     file_path = os.path.join(base_path, sample_dir, f"{sample_dir}_output", f"{sample_dir}_classification", f"{sample_dir}_annotated_cycles_files", "*.txt")
#     # Use glob to match all .txt files in the directory.
#     txt_files = glob.glob(file_path)

#     for txt_file in txt_files:
#         file_name = os.path.basename(txt_file).split('.')[0]
#         print("Processing file:", file_name)
        
#         with open(txt_file, 'r') as file:
#             for line in file:
#                 if line.startswith('Cycle='):
#                     # Remove 'Cycle=' and split the line by ';' then by '=' to create key-value pairs.
#                     cycle_data = dict(item.split('=') for item in line.split(';'))
#                     # Add the file name as the first entry in the dictionary.
#                     cycle_data['file_name'] = file_name
#                     # Append the dictionary to the list.
#                     all_cycle_data.append(cycle_data)

# df_cycles = pd.DataFrame(all_cycle_data)

# df_cycles['aa_barcode'] = df_cycles['file_name'].str.split('_', n=3).str[0]
# df_cycles['amplicon_number'] = df_cycles['file_name'].str.split('_', n=3).str[1]
# df_cycles = df_cycles[['file_name', 'aa_barcode', 'amplicon_number', 'Cycle', 'Copy_count', 'Length', 'IsCyclicPath', 'CycleClass', 'Segments' ]]

all_df_cycles = []
for sample_dir in sample_dirs:
    print("Processing", sample_dir)
    # Construct the file path for the amplicon classification profiles TSV file
    file_path = os.path.join(base_path, sample_dir, f"{sample_dir}_output", f"{sample_dir}_classification", f"{sample_dir}_annotated_cycles_files", "*.txt")
    # Use glob to match all .txt files in the directory.
    txt_files = glob.glob(file_path)
    
    for txt_file in txt_files:
        file_name = os.path.basename(txt_file).split('.')[0]
        print("Processing file:", file_name)
        all_cycle_data = []

        with open(txt_file, 'r') as file:
            for line in file:
                if line.startswith('Cycle='):
                    cycle_data = dict(item.split('=') for item in line.split(';'))
                    all_cycle_data.append(cycle_data)
        df_cycles = pd.DataFrame(all_cycle_data)
        df_cycles['file_name'] = file_name
        df_cycles['aa_barcode'] = df_cycles['file_name'].str.split('_', n=3).str[0]
        df_cycles['amplicon_number'] = df_cycles['file_name'].str.split('_', n=3).str[1]
        
        # Initialize DataFrame with required columns
        df_cycles = df_cycles.reindex(columns=['file_name', 'aa_barcode', 'amplicon_number', 'Cycle', 'Copy_count', 'Length', 'IsCyclicPath', 'CycleClass', 'Segments'])
        
        # Read segment data and process
        with open(txt_file, 'r') as file:
            contents = file.read()
        segment_lines = [line for line in contents.strip().split('\n') if line.startswith('Segment')]
        df_segments = pd.DataFrame([re.split(r'\t+', line) for line in segment_lines],
                                   columns=['Segment', 'Segment_Number', 'Chr', 'Start', 'End'])
        df_segments['Intervals'] = df_segments['Chr'] + ':' + df_segments['Start'] + '-' + df_segments['End']

        for index, row in df_cycles.iterrows():
            segments = row['Segments']
            numbers = re.findall(r'\d+', segments)  # extract numbers from the segment field
            df_segments_sub = df_segments[df_segments['Segment_Number'].isin(numbers)]
            df_segments_sub['Sort_Order'] = df_segments_sub['Segment_Number'].apply(lambda x: numbers.index(x))
            df_segments_sub = df_segments_sub.sort_values('Sort_Order')

            unique_intervals = df_segments_sub['Intervals'].drop_duplicates().tolist()
            cycle_interval = ','.join(unique_intervals)

            df_cycles.at[index, 'Intervals'] = cycle_interval
            # print(df_cycles)
        
        all_df_cycles.append(df_cycles)

big_df_cycles = pd.concat(all_df_cycles, ignore_index=True)

big_df_cycles['amplicon_index'] = big_df_cycles['file_name'].apply(lambda x: x.split('_')[-3])
big_df_cycles['aa_barcode'] = big_df_cycles['file_name'].apply(lambda x: '_'.join(x.split('_')[:-3]))
big_df_cycles['cycle_index'] = 'cycle' + big_df_cycles['Cycle'].astype(str)
big_df_cycles['cycle_barcode'] = big_df_cycles['aa_barcode'] + '_' + big_df_cycles['amplicon_index'] + '_' + big_df_cycles['cycle_index']

# Select columns after adding 'Intervals' column
big_df_cycles = big_df_cycles[['file_name', 'aa_barcode', 'amplicon_index', 'cycle_index', 'cycle_barcode', 'Copy_count','Length', 'IsCyclicPath', 'CycleClass', 'Segments', 'Intervals']]

# print(big_df_cycles)

final_string_7 = f"{formatted_string}_{manipulate_date}_{process}.csv"
print("save file path: ", save_path+final_string_7)
big_df_cycles.to_csv(save_path + final_string_7, index=False)
print(f"{process} file created successfully.")


log_file.close()
