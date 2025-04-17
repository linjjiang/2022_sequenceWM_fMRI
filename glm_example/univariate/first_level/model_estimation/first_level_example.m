%-----------------------------------------------------------------------
% SMGS first-level analysis example
% By: Linjing Jiang
% Date: 03/29/2023

% To do the first-level analysis, make sure you have the
% preproc_bold_mgs1.nii.gz and preproc_bold_mgs2.nii.gz files under the
% 'data' folder. You can get these files from Linjing directly (email me!)

clear
close all
clc

%% First define some parameters
rootDir = pwd; % root directory is the current directory
runs = 2; % number of runs
numScans = [348 348]; % number of scans -  make sure that this number is correct
TR = 1; % TR

%% Specify the model (similar to what we did when we check the design matrix)

% 0. set up directories
% data directory
dataDir = dir([rootDir,'/../data/']);
dataDir = [dataDir(1).folder '/'];
% output directory
outputDir = [rootDir '/output/']; % spm mat output directory
if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

% 1. Job specs
% Begin creating jobs structure
jobs{1}.spm.stats.fmri_spec.dir = cellstr(outputDir);
jobs{1}.spm.stats.fmri_spec.timing.units = "secs";
jobs{1}.spm.stats.fmri_spec.timing.RT = TR;
jobs{1}.spm.stats.fmri_spec.timing.fmri_t = 16; % You may want to change this
jobs{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;

% 2. Job sess, across runs
for rr = 1:runs
    % file name of the functional images
    funcs = ['preproc_bold_mgs',num2str(rr),'.nii'];
    
    % file name of the design matrix
    dms = ['design_matrix_mgs',num2str(rr),'.mat'];
    
    % file name of motion regressors
    mregs = ['motion_regressors_mgs',num2str(rr),'.txt'];
    
    % Select scans across runs
    % understand spm_select: http://jpeelle.net/mri/misc/spm_select.html
    if isempty(dir([dataDir,funcs]))
        % first unzip .nii.gz files
        gunzip([dataDir,funcs,'.gz'],dataDir)
    end
    % select scans
    files = cellstr(spm_select('ExtFPList',dataDir,['^' funcs], 1:2000)); 
    % check if the number of scans you provided is actuallly correct
    if length(files) ~= numScans(rr)
        error('Incorrect number of scans!')
    end
    % assign the functional data to the job structure
    jobs{1}.spm.stats.fmri_spec.sess(rr).scans = files;
    
    % Experimental conditions
    jobs{1}.spm.stats.fmri_spec.sess(rr).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    jobs{1}.spm.stats.fmri_spec.sess(rr).multi = {[dataDir,dms]};
    jobs{1}.spm.stats.fmri_spec.sess(rr).regress = struct('name', {}, 'val', {});
    
    % motion regressors
    jobs{1}.spm.stats.fmri_spec.sess(rr).multi_reg = {[dataDir,mregs]};
    
    % hpf
    jobs{1}.spm.stats.fmri_spec.sess(rr).hpf = 128;
end

% the rest of the job fields 
jobs{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
jobs{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
jobs{1}.spm.stats.fmri_spec.volt = 1;
jobs{1}.spm.stats.fmri_spec.global = 'None';
jobs{1}.spm.stats.fmri_spec.mthresh = 0.8;
jobs{1}.spm.stats.fmri_spec.mask = {''};
jobs{1}.spm.stats.fmri_spec.cvi = 'none';

% 3. model estimation options
jobs{2}.spm.stats.fmri_est.spmmat(1) = {[outputDir 'SPM.mat']};
jobs{2}.spm.stats.fmri_est.write_residuals = 0;
jobs{2}.spm.stats.fmri_est.method.Classical = 1;

% 4. define contrast: you want to specify your own contrast. Here I have
% three contrasts: stimulus activation, delay activation and resopnse
% activations
jobs{3}.spm.stats.con.spmmat(1) = {[outputDir 'SPM.mat']};

jobs{3}.spm.stats.con.consess{1}.tcon.name = 'stim';
jobs{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0];
jobs{3}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
jobs{3}.spm.stats.con.consess{2}.tcon.name = 'delay';
jobs{3}.spm.stats.con.consess{2}.tcon.weights = [0 1 0];
jobs{3}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
jobs{3}.spm.stats.con.consess{3}.tcon.name = 'resp';
jobs{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
jobs{3}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';
jobs{3}.spm.stats.con.delete = 0;

% navigate to output directory, estimate GLM
cd(outputDir)
spm_jobman('run',jobs)
