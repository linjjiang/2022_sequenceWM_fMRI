% Check design matrix in spm 
% By: Linjing, 10/26/22

% This is an important step before you run your first-level model. 
% Check if your design matrix looks correct!

% Before you run this script, make sure that you have the following:
% 1) Your design matrix ('.mat' file)
% 2) Number of scans in total
% 3) Motion regressors ('.txt' file)

% For this tutorial, I have already put all relevant files inside the
% 'data' folder, including:
% design matrix for run 1 of mgs: design_matrix_mgs1.mat
% design matrix for run 2 of mgs: design_matrix_mgs2.mat
% 6 motion regressors: motion_regressors_mgs1.txt,
% motion_regressors_mgs2.txt

% The files you don't need:
% preprocessed bold image: preproc_bold_mgs1.nii.gz, preproc_bold_mgs2.nii.gz
% You don't need actual brain data to check the design matrix! It means
% that you can do this step during the experimental design phase when the
% data has not yet been collected. If you do this, then you don't need
% motion regressors and you will need to change 'motion_or_not' value from
% 1 to 0.

% Alternatively, you can use spm batch and choose 'fmri_model_specification
% (design only) instead of using this script.

clear
close all
clc

% remember to add spm to your path if you are using spm from a local
% computer (change the path accordingly)
addpath(genpath('/usr/local/spm12'))
% If you are using spm on a HPC cluster, then you can do something like:
% module load spm

% Now add the data folder to the path
addpath('../data/')

%% Parameters to be defined

% scan parameter
rt = 1; % RT

% number of runs
runs = 2; 

% number of scans per run
nscan = [348 348];

% name of the multiple condition files
dmfilename = 'design_matrix_mgs*';

% name of the output directory: where you want to put the spm.mat output
outdirname = 'output/';

% whether to put into motion regressor
motion_or_not = 1; % add motion regressors
mtfilename = 'motion_regressors_mgs*.txt';
% Alternatively, assign 0 to motion_or_not if you don't have motion
% regressors yet.

%% Set up directories and file paths

% output directory
outdir = [pwd '/' outdirname]; % It will create an output directory within the current folder
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% path of the multiple-condition design matrix
dmfiles = dir(['../data/',dmfilename]);

% path of the motion regressors
if motion_or_not
mrfiles = dir(['../data/',mtfilename]);
end

%% Construct the GLM
job{1}.spm.stats.fmri_design.dir = {outdir};
job{1}.spm.stats.fmri_design.timing.units = 'secs';
job{1}.spm.stats.fmri_design.timing.RT = rt;
job{1}.spm.stats.fmri_design.timing.fmri_t = 16;
job{1}.spm.stats.fmri_design.timing.fmri_t0 = 8;
for rr = 1:runs
job{1}.spm.stats.fmri_design.sess(rr).nscan = nscan(rr);
job{1}.spm.stats.fmri_design.sess(rr).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
job{1}.spm.stats.fmri_design.sess(rr).multi = {[dmfiles(rr).folder '/' dmfiles(rr).name]};
job{1}.spm.stats.fmri_design.sess(rr).regress = struct('name', {}, 'val', {});
if motion_or_not
job{1}.spm.stats.fmri_design.sess(rr).multi_reg = {[mrfiles(rr).folder '/' mrfiles(rr).name]};
else
    job{1}.spm.stats.fmri_design.sess(rr).multi_reg = {''};
end
job{1}.spm.stats.fmri_design.sess(rr).hpf = 128;
end
job{1}.spm.stats.fmri_design.fact = struct('name', {}, 'levels', {});
job{1}.spm.stats.fmri_design.bases.hrf.derivs = [0 0];
job{1}.spm.stats.fmri_design.volt = 1;
job{1}.spm.stats.fmri_design.global = 'None';
job{1}.spm.stats.fmri_design.mthresh = 0.8;
job{1}.spm.stats.fmri_design.cvi = 'none';

%% Run the model
cd(outdir)
spm('defaults', 'FMRI');
spm_jobman('run', job);

%% Check the output in spm

% Now you should see an 'output' folder and a 'spm.mat' file under that
% folder. 
% Open up spm and use 'Review' function to check the output!
spm fmri