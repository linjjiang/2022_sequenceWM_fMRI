%-----------------------------------------------------------------------
% Job saved on 13-Feb-2023 14:22:27 by cfg_util (rev $Rev: 6942 $)
% spm SPM - SPM12 (7219)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
clear 
clc

%% Specify directories
rootDir = '/gpfs/projects/LeungGroup/vswmda/';
subjects = {'S5463_HL','S5475_HL','S5477_HL','S5499_HL','S5514_HL','S5515_HL'}; 
tasks = {'mgs','smgs'};
dm_type = {'run'};
con_type = {'stim','delay','response'};
dd = 1;
for tt = 1:2
    for cc = 1:3

spmDir = ['spm/',tasks{tt},'/',dm_type{1},'/'];
outputDir = [rootDir '2ndLevel/' tasks{tt} '/' dm_type{dd} '/' con_type{cc} '/']; % spm mat output directory

for ss = 1:length(subjects)
dataDir = [rootDir,subjects{ss},'/',spmDir];
files{ss,1} = spm_select('FPList',dataDir,['^con_000',num2str(cc),'.nii']);
end

%% Model specification
% Model specification: factorial design
matlabbatch{1}.spm.stats.factorial_design.dir = {outputDir};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = files;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%% Model estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = {[outputDir 'SPM.mat']};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

%% Contrast manager
matlabbatch{3}.spm.stats.con.spmmat = {[outputDir 'SPM.mat']};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = con_type{cc};
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

%% Estimation
% navigate to output directory, specify and estimate GLM
cd(outputDir)
spm_jobman('run',matlabbatch)

    end
end
