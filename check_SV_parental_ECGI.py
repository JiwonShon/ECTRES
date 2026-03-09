# /data/ECTRES/results/aaSuite_germline_ms/v1.3.8/GRCh37/minCN4.5/cnsizeMin50000/-1X/calls/ECTRES-ECGI1-0001-TPX-A01-WGS-1ST985/ECTRES-ECGI1-0001-TPX-A01-WGS-1ST985_output/ECTRES-ECGI1-0001-TPX-A01-WGS-1ST985_classification/ECTRES-ECGI1-0001-TPX-A01-WGS-1ST985_SV_summaries/ECTRES-ECGI1-0001-TPX-A01-WGS-1ST985_amplicon7_SV_summary.tsv

# singularity shell -B /mnt/NAS3/home/jiwon/HL-NF/scratch/ECTRES:/data/ECTRES/ aasuite_v1.3.8.sif
 # singularity shell -B /mnt/NAS3/home/jiwon/p:/data/ECTRES/ aasuite_v1.3.8.sif

#python3 check_SV_support.py -i /mnt/NAS3/home/jiwon/HL-NF/scratch/ECTRES/results/aaSuite_germline_ms/v1.3.8/GRCh37/minCN4.5/cnsizeMin50000/-1X/calls/ECTRES-ECGI1-0001-TPX-A01-WGS-1ST985/ECTRES-ECGI1-0001-TPX-A01-WGS-1ST985_output/ECTRES-ECGI1-0001-TPX-A01-WGS-1ST985_classification/ECTRES-ECGI1-0001-TPX-A01-WGS-1ST985_SV_summaries/ECTRES-ECGI1-0001-TPX-A01-WGS-1ST985_amplicon7_SV_summary.tsv --bam /mnt/NAS3/home/jiwon/HL-NF/scratch/ECTRES/results/align_dna/GRCh37/ECTRES-ECGI1-0001-TPX-A23-WGS-3RT414/applyBqsr/ECTRES-ECGI1-0001-TPX-A23-WGS-3RT414.realn.mdup.bqsr.bam --ref GRCh37 --tolerance 200 --window 500
