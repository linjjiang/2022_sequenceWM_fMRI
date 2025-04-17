% Prepare beta maps for MVPA
% By: Linjing Jiang
% Time: 03/26/2024

% This script is to prepare beta maps for MVPA or RSA

% Input:

% beta_*.nii: first-level GLM results from SPM in specific subject and
% session folder

% beta_label.mat: output from the 'get_beta_labels.m' that indicates which
% beta maps are relevant, the run number and condition for the beta maps.
% This mat file should be saved in the same folder as in beta*.nii files
% for now, beta maps and labels are saved in the following folder, e.g.,
% for subject f16 and session S5824
% '/mnt/Data1/linjdata1/vswmda/scan_data/spm/output/preproc_st_spatiotemporal/f16/S5824/smgs'

% roi: for now (03/26/2024), ROIs are saved as .nii files inside the
% '/mnt/Data1/linjdata1/vswmda/atlases/final_atlas_roi' folder
% for subject f16, the ROIs are the following
% mgs_*(roi name)_f16.nii


% Outputs:

% output folder is:
% '/mnt/Data1/linjdata1/vswmda/scan_data/rsa/beta_run_rsa/f16/beta' for now for
% subject f16

% we will output extracted beta maps for each ROI, and name it as
% beta_session_run_condition_roi.nii
% for example, for f16 session S5824, if we extract beta map for delay_3WL
% condition and for roi v1, it will be named as following:
% beta_S5824_R1_delay_3WL_v1.mat

% Specifically, we will define a function, extract_roi_data, to extract
% beta maps using specific ROIs. This function is saved in the following
% path: /mnt/Data1/linjdata1/vswmda/script/make_roi

function prep_beta_map_for_mvpa_mgsROI(rootDir,rootOutDir,roiDir,roi_pfix,roi_names,idx)
%% clean up the environment and set up the directories
addpath(genpath('/mnt/Data1/linjdata1/vswmda/script/'))
addpath(genpath('/usr/local/spm12'))

% clear
% close all
% clc

% % root data directory (change as needed)
% rootDir = '/mnt/Data1/linjdata1/vswmda/scan_data/spm/output/preproc_st_spatiotemporal_censor';
% 
% % root output directory
% rootOutDir = '/mnt/Data1/linjdata1/vswmda/scan_data/rsa/beta_run_rsa_censor/';

% subject directory. change if necessary.
subDir = dir(fullfile(rootDir,'f*'));

% task
task_id = 'smgs'; % task name

% % roi directory (change as needed)
% roiDir = '/mnt/Data1/linjdata1/vswmda/atlases/final_atlas_roi';

%%
for sb = idx%1:length(subDir) % for each subject

    sub_id = subDir(sb).name; % subject id
    
        % load ROIs
        for rr = 1:length(roi_names)
   roi_files(rr) = dir(fullfile(roiDir,[roi_pfix,roi_names{rr},'_',sub_id,'.nii'])); % all the roi files %'mgs_*_'
        end
        
% % get roi names
% roi_names = struct2cell(roi_files); roi_names = roi_names(1,:);
% roi_names = cellfun(@(C) erase(C,{roi_pfix,sub_id,'.nii'}),roi_names,'UniformOutput',false);

    sessDir = dir(fullfile(subDir(sb).folder,subDir(sb).name,'S*')); % session directories
    
    % output directory - we will have one output folder for each
    % participant (combining multiple sessions
    outputDir = fullfile(rootOutDir,sub_id,'/beta'); % extracted beta maps output directory
    
    if ~exist(outputDir,'dir') % if the output directory is empty, we make it
        mkdir(outputDir)
    end
    
    for ss = 1:length(sessDir) % for each session under that subject
        
        sess_id = sessDir(ss).name; % session id
        
        clearvars -except idx rootDir rootOutDir subDir sessDir outputDir sb ss sub_id sess_id task_id roiDir roi_files roi_names roi_pfix
        
        % get the actual data directory
        dataDir = fullfile(rootDir,sub_id,sess_id,task_id); % spm mat output directory
        if ~isempty(dir(fullfile(dataDir,'beta_label.mat'))) % if the data directory is not empty (there is already a SPM result)
            
            % we first load the relevant beta maps
            load(fullfile(dataDir,'beta_label.mat'))
            
            % concatenate beta maps of all conditions and all task periods
            beta_tbl = [stim_beta;delay_beta;response_beta];
            
            % extract data from each beta map using specific roi
            for ii = 1:length(roi_files) % for each roi
                % define roi path
                roi_path = fullfile(roi_files(ii).folder,roi_files(ii).name);
                
                % get roi name
                roi_name = roi_names{ii};
                
                % then we want to extract specific beta maps from each
                % condition and run
                for bb = 1:size(beta_tbl,1) % for each unique condition and run
                    
                    % we get the corresponding beta map names
                    beta_name = beta_tbl.beta_name(bb);
                    beta_path = cellfun(@(C) fullfile(dataDir,C),beta_name,'UniformOutput',false);
                    
                    % extract beta maps data using roi
                    roi_beta = extract_roi_data(roi_path,beta_path);
                    
                    % get condition name
                    cond_name = beta_tbl.condition_name{bb};
                    
                    % get run name
                    run_name = ['R' num2str(beta_tbl.run_index(bb))];
                    
                    % new file name for extracted beta maps
                    beta_name_new = ['beta_' sess_id '_' run_name '_' cond_name '_' roi_name '.mat'];
                    
                    % save the new beta data
                    save(fullfile(outputDir,beta_name_new),'roi_beta')
                    clear roi_beta
                end % for each beta map
            end % for each roi        
        end % if the data directory is not empty (has beta maps)
    end % for each session/data directory
end % for each participant

end