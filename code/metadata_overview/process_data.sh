#!/bin/bash
#SBATCH --job-name=dada2
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p cu
#SBATCH --mem=30G
#SBATCH -o output/dada2_%a.out
#SBATCH -e output/dada2_%a.err
#SBATCH --array=2,7,8
#SBATCH --exclude=cu01

module load R/4.1.3

# 安全读取study_id和路径        
study_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" study_id.txt | tr -d '\r\n')
file_path="../../data/projects/${study_id}"
path=/home/dongbiao/software/parallel-20250222/bin

# makedirs
# mkdir ${file_path}/intermediate ${file_path}/rawdata ${file_path}/results ${file_path}/temp

# download rawdata
# prefetch --option-file ${file_path}/sample.txt --output-directory ${file_path}/sra
# ls ${file_path}/sra/*/*.sra | awk -F'/' '{print $8}' | awk -F'.' '{print $1}' > ${file_path}/sample.txt

threads=8
# 并行运行 fastq-dump
# cat "${file_path}/sample.txt" | ${path}/parallel -j ${threads} --bar \
    # "fastq-dump --split-3 ${file_path}/sra/{}/{}.sra --outdir ${file_path}/rawdata"

# Run DADA2
ls ${file_path}/rawdata/ > ${file_path}/run_sample.txt

Rscript process_forwards.R ${file_path} trim
