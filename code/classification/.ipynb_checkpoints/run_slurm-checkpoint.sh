#!/bin/bash
#SBATCH --job-name=xgboost
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p cu
#SBATCH --mem=250G
#SBATCH -o xgboost.out
#SBATCH -e xgboost.err
#SBATCH --exclude=cu01

module load miniconda/4.9.2
source activate jupyter_notebook

python Hyperparameters_Xgboost.py