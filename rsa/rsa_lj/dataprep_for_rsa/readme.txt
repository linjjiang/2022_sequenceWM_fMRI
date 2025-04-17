Prepare beta maps for running RSA
Date: 09/11/2024
Contact: Linjing Jiang, linjjiang07@gmail.com

- /mnt/Data1/linjdata1/vswmda/script/used_for_final_analysis/rsa_lj/dataprep_for_rsa
	prepare data for rsa using full GLM (all 8 spatiotemporal conditions as different regressors)
	In MATLAB, extract beta maps using ROIs and prepare data structure for running RSA in python.
	Using the following types of ROIs
- using thresholded MGS ROI, made by overlapping atlas ROI and conjunctive activation of mgs1 and mgs2 (two sessions)
	- all events (stimulus + delay + response), p < 0.05 voxel-wise, cluster size > 50
- using atlas ROI: 
	- v1/v2: wang probability atlas 25%
	- lateral parietal: wang probability atlas 15% (why 15% -- we want to be more liberal as we will overlap this with mgs activation)
	- fef: wang probability atlas 25%
	- area4: julich atlas probability 50%
	- sfg, mfg, ifg: from MD network
	
1) main_scripts: main scripts for data prep. You may need to change the path if needed. Including:
- main_rsa_prep_atlasROI.m : prepare data using atlas ROI
- main_rsa_prep_mgsROI_thresMGS.m : prepare data using thresholded MGS ROI
- main_rsa_prep_mgsROI_topN.m : prepare data using topN MGS ROI
- main_rsa_prep_mgsROI_left_topN.m AND matlab_left_topN.slurm : prepare data using topN ROI, left hemisphere, dilated/closed
- main_rsa_prep_mgsROI_right_topN.m AND matlab_right_topN.slurm : prepare data using topN ROI, right hemisphere, dilated/closed

The rest of the scripts are not used for final analysis and stored in not_used folder:
- main_rsa_prep_mgsROI_left_thresMGS.m AND matlab_left_thresMGS.slurm : prepare data using thresholded MGS ROI, left hemisphere, not dilated/closed
- main_rsa_prep_mgsROI_right_thresMGS.m AND matlab_right_thresMGS.slurm : prepare data using thresholded MGS ROI, right hemisphere, not dilated/closed
- main_rsa_prep_mgsROI_left_thresMGS_close.m AND matlab_left_thresMGS_close.slurm : prepare data using thresholded MGS ROI, left hemisphere, dilated/closed
- main_rsa_prep_mgsROI_right_thresMGS_close.m AND matlab_right_thresMGS_close.slurm : prepare data using thresholded MGS ROI, right hemisphere, dilated/closed


2) functions: this folder contains the following functions in order:
get_beta_labels.m: get labels of the beta map for exp1
get_beta_labels_exp2.m: get labels of the beta map for exp2
prep_beta_map_for_mvpa_atlasROI: prepare beta maps for atlas ROI
prep_beta_map_for_mvpa_mgsROI: prepare beta maps for mgs ROI including thresholded MGS and top N MGS
rsa_data_prep_beta_by_run: constuct data structure that is input into python script for doing RSA for exp 1
rsa_data_prep_beta_by_run_exp2: constuct data structure that is input into python script for doing RSA for exp2
backup_functions: just in case the previous functions do not work,  you can check the scripts in the backup_functions folders and debug

3) old: old scripts, may be deleted

4) slurm_logs: logs of the slurm jobs submitted to HPC


Step 4: Run RSA
/mnt/Data1/linjdata1/vswmda/script/used_for_final_analysis/rsa_lj/xxxx:
1. rsa_atlas_ROI_updated: RSA using atlas ROI
2. rsa_mgs1_conj_mgs2_0.05_50: RSA using thresholded MGS ROI
3. rsa_mgs1_conj_mgs2_topN_updated: RSA using top N ROI
4. rsa_mgs1_conj_mgs2_topN_left/right: RSA using top N ROI left and right hemispheres

