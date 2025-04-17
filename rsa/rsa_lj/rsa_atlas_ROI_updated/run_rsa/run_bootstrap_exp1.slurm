#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=96
#SBATCH --time=01:00:00
#SBATCH --partition=long-96core
#SBATCH --job-name=exp1
#SBATCH --output=run_bootstrap_exp1_output.txt

module load anaconda/3
source activate glmenv

python run_bootstrap_exp1.py
