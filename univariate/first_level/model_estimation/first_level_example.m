%-----------------------------------------------------------------------
% SMGS first-level analysis
clear
clc

%% First define some parameters
rootDir = '/gpfs/projects/LeungGroup/vswmda/';
subjects = {'S5463_HL'};%{'S5475_HL','S5477_HL','S5499_HL','S5514_HL','S5515_HL'}; 

tasks = {'mgs'};
dm_type = 'run';

if strcmp(tasks{1},'smgs')
% number of scans and TR
runs = 4; % number of runs
numScans = [419 419 419 419];
disScans = [0 0 0 0]; % number of scans discard during preprocessing
elseif strcmp(tasks{1},'mgs')
runs = 2; % number of runs
numScans = [348 348];
disScans = [0 0]; % number of scans discard during preprocessing
end

numScans = numScans - disScans;
TR = 0.99;

ESTIMATE_GLM = 1; % Whether to estimate GLM

%% For each subject
for sub = subjects
    for tsk = tasks
        % directories
        dataDir = dir([rootDir,'/derivatives/halfpipe/sub-',sub{1}(1:5),'HL/func/']);
        dataDir = [dataDir(1).folder '/']; % data directory (nii)
        spmDir = ['spm/',tsk{1},'/',dm_type,'/'];
        outputDir = [rootDir sub{1} '/' spmDir]; % spm mat output directory

        % 1. Job specs
        % Begin creating jobs structure
        jobs{1}.spm.stats.fmri_spec.dir = cellstr(outputDir);
        jobs{1}.spm.stats.fmri_spec.timing.units = "secs";
        jobs{1}.spm.stats.fmri_spec.timing.RT = TR;
        jobs{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        jobs{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;

        % 2. Job sess, across runs
        for rr = 1:runs
            % functional images
            funcs = ['sub-',sub{1}(1:5),'HL_task-',tsk{1},'_run-',num2str(rr),'_setting-preproc_bold.nii'];

            % design matrix
            dmDir = ['dm_',dm_type,'_',tsk{1},num2str(rr),'.mat'];

            % motion regressors
            mregs = ['sub-',sub{1}(1:5),'HL_task-',tsk{1},'_run-',num2str(rr),'_setting-preproc_desc-confounds_regressors.txt'];

            % Select scans across runs
            % understand spm_select: http://jpeelle.net/mri/misc/spm_select.html
            if isempty(dir([dataDir,funcs]))
                % unzip .nii.gz files
                gunzip([dataDir,funcs,'.gz'],dataDir)
            end
            % select scans
            files = cellstr(spm_select('ExtFPList',dataDir,['^' funcs], 1:500)); %numScans(rr)
%             % check if the number of scans equal the number we expect
%             if length(files) ~= numScans(rr)
%                 error('Incorrect number of scans!!')
%             end

            jobs{1}.spm.stats.fmri_spec.sess(rr).scans = files;

            % Conditions
            jobs{1}.spm.stats.fmri_spec.sess(rr).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
            jobs{1}.spm.stats.fmri_spec.sess(rr).multi = {[rootDir sub{1} '/' dmDir]};
            jobs{1}.spm.stats.fmri_spec.sess(rr).regress = struct('name', {}, 'val', {});

            % motion regressors
            jobs{1}.spm.stats.fmri_spec.sess(rr).multi_reg = {[dataDir mregs]};

            % hpf
            jobs{1}.spm.stats.fmri_spec.sess(rr).hpf = 128;
        end

        % the rest of the job fields (what are those fields?)
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

        % 4. define contrast
        jobs{3}.spm.stats.con.spmmat(1) = {[outputDir 'SPM.mat']};
        if strcmp(dm_type,'run')
            jobs{3}.spm.stats.con.consess{1}.tcon.name = 'stim';
            jobs{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0];
            jobs{3}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
            jobs{3}.spm.stats.con.consess{2}.tcon.name = 'delay';
            jobs{3}.spm.stats.con.consess{2}.tcon.weights = [0 1 0];
            jobs{3}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
            jobs{3}.spm.stats.con.consess{3}.tcon.name = 'resp';
            jobs{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
            jobs{3}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';
        elseif strcmp(dm_type,'order')
            jobs{3}.spm.stats.con.consess{1}.tcon.name = 'ord1';
            jobs{3}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0];
            jobs{3}.spm.stats.con.consess{1}.tcon.sessrep = 'replsc';
            jobs{3}.spm.stats.con.consess{2}.tcon.name = 'ord2';
            jobs{3}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 0];
            jobs{3}.spm.stats.con.consess{2}.tcon.sessrep = 'replsc';
            jobs{3}.spm.stats.con.consess{3}.tcon.name = 'ord3';
            jobs{3}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 0];
            jobs{3}.spm.stats.con.consess{3}.tcon.sessrep = 'replsc';
            jobs{3}.spm.stats.con.consess{4}.tcon.name = 'ord4';
            jobs{3}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1];
            jobs{3}.spm.stats.con.consess{4}.tcon.sessrep = 'replsc';
            jobs{3}.spm.stats.con.consess{5}.tcon.name = 'ordleft';
            jobs{3}.spm.stats.con.consess{5}.tcon.weights = [0 1 1 0];
            jobs{3}.spm.stats.con.consess{5}.tcon.sessrep = 'replsc';
            jobs{3}.spm.stats.con.consess{6}.tcon.name = 'ordright';
            jobs{3}.spm.stats.con.consess{6}.tcon.weights = [1 0 0 1];
            jobs{3}.spm.stats.con.consess{6}.tcon.sessrep = 'replsc';
            jobs{3}.spm.stats.con.consess{7}.tcon.name = 'ordleft-ordright';
            jobs{3}.spm.stats.con.consess{7}.tcon.weights = [-1 1 1 -1];
            jobs{3}.spm.stats.con.consess{7}.tcon.sessrep = 'replsc';
        end
        jobs{3}.spm.stats.con.delete = 0;

        % navigate to output directory, specify and estimate GLM
        cd(outputDir)
        spm_jobman('run',jobs)


        %         if ESTIMATE_GLM == 1
        %             load SPM;
        %             spm_spm(SPM);
        %             %             matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        %             %             matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        %             %             matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        %         end
    end
end