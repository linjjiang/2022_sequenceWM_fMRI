#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=96
#SBATCH --time=04:00:00
#SBATCH --partition=short-96core
#SBATCH --job-name=exp2
#SBATCH --output=run_bootstrap_exp2_output.txt

module load anaconda/3
source activate glmenv

python run_bootstrap_exp2.py
