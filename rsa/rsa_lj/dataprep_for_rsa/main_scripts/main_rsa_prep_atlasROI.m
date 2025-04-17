clear all
close all
clc

cd('/mnt/Data1/linjdata1/vswmda/script/rsa_lj/rsa_use_full_GLM/prep_data_atlasROI')

% spm output directory (data directory)
rootDir = '/mnt/Data1/linjdata1/vswmda/scan_data/spm/output/preproc_st_spatiotemporal_censor'; 

% rsa output directory
rootOutDir = '/mnt/Data1/linjdata1/vswmda/scan_data/rsa/full_GLM_atlas_roi/';

% roi directory (change as needed)
roiDir = '/mnt/Data1/linjdata1/vswmda/atlases/final_atlas_roi/rsa_atlas_roi'; 
% roi_pfix = 'mgs_'; % 'mgs_dil_'; 'mgs_75_'; 'mgs_55_dil_';
% roi_names = {'v1', 'v2', 'area4', 'fef', 'sfg', 'mfg', 'ifg', 'ips',...
%     'ips0', 'ips1', 'ips2', 'ips3', 'ips4', 'ips5', 'spl1'};

subDir = dir(fullfile(rootDir,'f*')) % subject directory. change if necessary.

% if this is already run for previous sets of analysis, you can comment
% this
get_beta_labels(rootDir,1:6)
get_beta_labels_exp2(rootDir,7:9)

prep_beta_map_for_mvpa_atlasROI(rootDir,rootOutDir,roiDir,1:9)

rsa_data_prep_beta_by_run(rootOutDir,1:6)
rsa_data_prep_beta_by_run_exp2(rootOutDir,7:9)


% %%
% clear all
% close all
% clc
% 
% cd('/mnt/Data1/linjdata1/vswmda/script/rsa_lj/complex_rsa_model')
% 
% % spm output directory (data directory)
% rootDir = '/mnt/Data1/linjdata1/vswmda/scan_data/spm/output/preproc_nsm_spatiotemporal_censor'; 
% % rsa output directory
% rootOutDir = '/mnt/Data1/linjdata1/vswmda/scan_data/rsa/beta_run_nsm_rsa_censor/';
% 
% get_beta_labels(rootDir)
% prep_beta_map_for_mvpa(rootDir,rootOutDir,roiDir,roi_pfix,roi_names)
% rsa_data_prep_beta_by_run(rootOutDir)