#!/bin/bash

#SBATCH --job-name=function_predict
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p cu
#SBATCH -w cu02
#SBATCH --mem=200G
#SBATCH -o out.out
#SBATCH -e out.err

module load miniconda/4.9.2
source activate picrust2
picrust2_pipeline.py -s /home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/projects/feces_seq_16S.fasta -i /home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/projects/table_3603_2339.biom -o/home/dongbiao/word_embedding_microbiome/programe_test/dietary_fber/data/picrust2 -p 28 --in_traits GO,EC,KO,PFAM,CAZY,BIGG --stratified --verbose
