#!/bin/bash
#SBATCH --job-name=right_topN       # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH -p gpu
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
##SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=04:00:00          # total run time limit (HH:MM:SS)
##SBATCH --mail-type=all          # send email on job start, end and fault
##SBATCH --mail-user=<YourNetID>@princeton.edu
#SBATCH --output right_topN.log

#module purge
module load matlab spm shared

cd /gpfs/projects/LeungGroup/hc_backup/vswmda/script/rsa_lj/rsa_use_full_GLM/prep_data_mgsROI
matlab -singleCompThread -nodisplay -nosplash -r main_rsa_prep_mgsROI_right_topN
