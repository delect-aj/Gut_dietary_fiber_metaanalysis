#!/bin/bash

#SBATCH --job-name=cluster_features
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p cu
#SBATCH --mem=30G
#SBATCH -o output/cluster_features_%a.out
#SBATCH -e output/cluster_features_%a.err
#SBATCH --array=1-9
#SBATCH --exclude=cu01

module load miniconda/4.9.2 
source activate qiime2-amplicon-2024.2

study_id=$(sed -n "${SLURM_ARRAY_TASK_ID}p" study_id.txt | tr -d '\r\n')
file_path="../../data/projects/${study_id}"

# ASVs_counts.tsv to ASVs_counts.biom
biom convert \
  -i ${file_path}/results/ASVs_counts.tsv \
  -o ${file_path}/results/ASVs_counts.biom \
  --table-type="OTU table" \
  --to-hdf5

qiime tools import \
  --input-path ${file_path}/results/ASVs.fa \
  --output-path ${file_path}/results/rep-seqs.qza \
  --type 'FeatureData[Sequence]'

qiime tools import \
  --input-path ${file_path}/results/ASVs_counts.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path ${file_path}/results/table.qza
  
qiime vsearch cluster-features-closed-reference \
    --i-sequences ${file_path}/results/rep-seqs.qza \
    --i-table ${file_path}/results/table.qza \
    --i-reference-sequences /beegfs/dongbiao/greengene2/2024.09.backbone.full-length.fna.qza \
    --p-perc-identity 0.97 --p-threads 12 \
    --p-strand both \
    --o-clustered-table ${file_path}/results/clustered-table \
    --o-clustered-sequences ${file_path}/results/clustered-sequences \
    --o-unmatched-sequences ${file_path}/results/unmatched-sequences

qiime tools export \
  --input-path ${file_path}/results/clustered-table.qza \
  --output-path ${file_path}/results/extracted-feature-table
  
