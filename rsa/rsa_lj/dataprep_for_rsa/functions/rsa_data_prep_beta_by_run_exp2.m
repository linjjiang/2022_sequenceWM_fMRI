% Prepare data for RSA using beta maps per run per condition
% By: Linjing Jiang
% Date: 03/26/2024

% In this script, I prepared data for RSA later
% The data was collected in 2023.7 - 2023.12.
% 5 participant has 3 sessions, one session consists of visuomotor tasks
% Another two sessions consist of memory tasks
% In total, we ran 12 runs of smgs tasks, each run has 16 trials, resulting
% in 192 trials in total.

% The main analysis focuses on decoding eight spatiotemporal patterns, each
% consisting of 24 trials.
% These eight patterns include:
% 3WL: 3dva, wide, L->R (0,0,0)
% 3TL: 3dva, tall, L->R (0,1,0)
% 3WR: 3dva, wide, R->L (0,0,1)
% 3TR: 3dva, tall, R->L (0,1,1)
% 5.5WL: 5.5dva, wide, L->R (1,0,0)
% 5.5TL: 5.5dva, tall, L->R (1,1,0)
% 5.5WR: 5.5dva, wide, R->L (1,0,1)
% 5.5TR: 5.5dva, tall, R->L (1,1,1)

% The input is:
% Beta map per experimental condition per run, that is originally
% produced by first-level GLM in SPM using the spatiotemporal patterns as
% regressors (preproc_st_spatiotemporal), and then extracted using the
% '/mnt/Data1/linjdata1/vswmda/script/rsa_lj/prep_beta_map_for_mvpa.m' and
% '/mnt/Data1/linjdata1/vswmda/script/rsa_lj/get_beta_labels.m'

% Those beta maps are in the following location:
% '/mnt/Data1/linjdata1/vswmda/scan_data/rsa/beta_run_rsa_censor/f16/beta'
% for subject f16

% Those beta maps have naming like 
% 'beta_S5847_R4_response_5.5TR_v2.mat'
% In this example, we have the extracted beta signals from session S5847,
% run 4 (not the actual run number, just the order of the runs in the GLM
% regressors), condition response, 5.5dva, tall, and right first, ROI v2


% The outputs are:
% 'responsePatterns': data ready for RSA
% This is a 1 x 1 structure, containing 
% responsePatterns.stim
% responsePatterns.delay
% responsePatterns.response
% corresponding to different task epochs(stimulus, delay, and response)
% output directory:  '/mnt/Data1/linjdata1/vswmda/scan_data/rsa/beta_run_rsa_censor/f16/'

% for each task epoch, e.g., responsePatterns.stim
%   it is a 1 x N structure, where N is the number of ROIs x number of
%   sessions x number of participants. For example, if there are 7 ROIs, 1
%   participants and 2 sessions per participant, then we get a 1 x 14
%   structure. 
%   Within each structure, there are several fields:
%   'name': ROI name | session | participant id
%   'data': a V x R x C matrix, where V is the number of voxels in that
%   ROI, R is the number of runs (not trials, because we are using run-wise 
%   beta maps derived from first-level GLM), C is the number of conditions. For 
%   example, if there are 400 voxels in that ROI, 3 runs of beta maps per 
%   condition, and 8 conditions, then we get a 400 x 3 x 8 matrix.

function rsa_data_prep_beta_by_run_exp2(rootDir,idx)
%% 1. First of all, specify which dataset(s) you want to run
% clear
% close all
% clc

addpath(genpath('/mnt/Data1/linjdata1/vswmda/script/'))

% % root directory
% rootDir = '/mnt/Data1/linjdata1/vswmda/scan_data/rsa/beta_run_rsa_censor/';

% subject directory. change if necessary.
subDir = dir(fullfile(rootDir,'f*'));

for sb = idx %1:length(subDir)

% specify subject id
sub_id = subDir(sb).name;

out_dir = fullfile(rootDir,sub_id); % output directory
% use a new set of ROI defined in March 2024 using Wang+AAL+Julich atlas and functional activation of mgs
% use run-wise beta maps masked with each ROI, instead of using
% preprocessed bold timecourse as before

data_dir = fullfile(rootDir,sub_id,'beta'); % data directory (where the extracted beta maps are stored)

% if ~isfolder(out_dir)
%  mkdir(out_dir)
% end

%% 2. We then load all the beta maps from the data directory
% make sure to clear those variables first if they are in the workspace
% already
clear responsePatterns*

% load all the beta maps from the data directory
beta_files = dir(fullfile(data_dir,'beta*'));

% get the file names
beta_names = struct2cell(beta_files); beta_names = beta_names(1,:);
beta_names_all = split(beta_names,"_"); beta_names_all = squeeze(beta_names_all);
beta_names_all(:,6) = cellfun(@(C) erase(C,'.mat'),beta_names_all(:,6),'UniformOutput',false); % erase '.mat' at the end to get pure roi label

% responsePatterns.delay, stim, resp: task epoch
% under each, we have 
% 'name': ROI name | participant id | session | run
% 'rawdata': a V x C matrix, where V is the number of voxels in that
    %  ROI, C is the number of conditions. For 
    %   example, if there are 400 voxels in that ROI, 
    %   and 8 conditions, then we get a 400 x 8 matrix. Each cell contains
    %   the raw data retrieved from the beta maps.
% 'data': we have removed nans from the rawdata.

% let's figure out the unique file names

% all session names
sess = beta_names_all(:,2);

% unique session names
uniq_sess = unique(sess);

% all run names
run = beta_names_all(:,3);

% unique run names
uniq_run = unique(run);

% all task period
epoch = beta_names_all(:,4);

% unique task period
uniq_epoch = unique(epoch);

% all condition
cond = beta_names_all(:,5);

% unique condition
uniq_cond = unique(cond);

% intended order of conditions
good_cond = {'CWL1','CTL1','CWR1','CTR1',...
        'CWL2','CTL2','CWR2','CTR2',...
        'NWL1','NTL1','NWR1','NTR1',...
        'NWL2','NTL2','NWR2','NTR2'};
% good_cond = {'CWL','CTL','CWR','CTR',...
%         'CWL','CTL','CWR','CTR'};
 
% all roi names
roi = beta_names_all(:,6);

% unique roi names
uniq_roi = unique(roi);

for ee = 1:length(uniq_epoch) % for each task epoch, we will get a different variable responsePatterns
    curr_epoch = uniq_epoch{ee}; % current task epoch
    
    kk = 1; % an counter for field names
for ii = 1:length(uniq_roi) % for each roi
    for ss = 1:length(uniq_sess) % for each session
        for rr = 1:length(uniq_run) % for each run
        
        % for each roi, subject, and session, we will get a unique field
        % under the structure responsePatterns
        curr_roi = uniq_roi{ii}; % current roi
        curr_sub = sub_id; % current subject
        curr_sess = uniq_sess{ss}; % current session
        curr_run = uniq_run{rr}; % current run
        
        % we need to get the run idx in the beta map matrix
        run_idx = str2num(erase(curr_run,'R'));
       
        % name of the field        
        fn = [curr_roi ' | ' curr_sub ' | ' curr_sess ' | ' curr_run];
        
        % get a random beta map for current roi
        a_random_beta_for_curr_roi = find(strcmp(roi,curr_roi)); 
        a_random_beta_for_curr_roi = a_random_beta_for_curr_roi(1);
        a_beta = load(fullfile(beta_files(a_random_beta_for_curr_roi).folder,beta_files(a_random_beta_for_curr_roi).name),'roi_beta');
        a_beta = a_beta.roi_beta;

        clear beta_VC
        % initialize beta_VC matrix
        beta_VC = nan(length(a_beta),length(uniq_cond));
        
        % for the specific roi, subject, session, and run, we will
        % retrieve the beta data as a V x C matrix (voxel by
        % condition)
        for cc = 1:length(uniq_cond) % for each condition

                clear curr_beta
                curr_cond = uniq_cond{cc}; % current conditionn
                 
                % get the condition idx
                % we already know the intended order of conditions, stored
                % in 'good_cond' variable. We just need to find out the
                % index of the current condition among the good condition
                % variable
                cond_idx = strcmp(good_cond,curr_cond);
                
                % get the index of the beta maps regarding 'beta_files'
                % variable
                idx = find(strcmp(epoch,curr_epoch) & ...
                           strcmp(roi,curr_roi) & ...
                           strcmp(sess,curr_sess) & ...
                           strcmp(cond,curr_cond) & ...
                           strcmp(run,curr_run)); % there should only be one index
                       
                % for some participants, not each session has all the conditions
                % and runs. For those, we will store them as nans.
       
                if ~isempty(idx)
                    if length(idx) == 1
                % load the specific beta map
                curr_beta = load(fullfile(beta_files(idx).folder,beta_files(idx).name),'roi_beta');
                curr_beta = curr_beta.roi_beta;
                    elseif length(idx) > 1
                        error('More than one beta map per run per condition per session per roi per subject!')
                    end
                
                % store the beta map. Be careful that we need to transpose
                % the data because curr_beta is a row vector
                beta_VC(:,cond_idx) = curr_beta';  
                
                % if current condition is spatial set 1 or 2
                if strcmp(curr_cond(end),'1') % if 3 dva
                    curr_set = 'P1';
                else
                    curr_set = 'P2';
                end
                
                end    
            end
% update the field name with eccentricity
        fn = [curr_roi ' | ' curr_sub ' | ' curr_sess ' | ' curr_run ' | ' curr_set];
        
%                 % we also want to get rid of columns with all NANS (which
%                 % means during that run, there is no such condition, e.g.,
%                 % eccentricity)
%                 % That is, as long as there is a run that contains data
%                 % from a specific condition, we will include that run's
%                 % data (no matter which run it is).
%                 % for example, normally there are three runs of 3 dva and
%                 % three runs of 5.5 dva stimuli. Then we will reduce the
%                 % data matrix to 3 columns from 6 columns.
%                  idx_has_data = squeeze(~all(ismissing(beta_VRC),1));
%                  % just in case the number of runs are not balanced across
%                  % conditions, i.e., some participants did not complete all
%                  % the runs in that session so some eccentricities have
%                  % less number of runs than others. Spatiotemporal
%                  % conditions are balanced though.
%                  beta_VRC_has_data = nan(size(beta_VRC,1),max(sum(idx_has_data,1)),size(beta_VRC,3)); % Number of voxels, maximum number of runs of data, number of conditions
% 
%                  for jj = 1:size(idx_has_data,2)
%                      row_idx = find(idx_has_data(:,jj)==1);
%                      beta_VRC_has_data(:,1:length(row_idx),jj) = beta_VRC(:,row_idx,jj);
%                  end
                                       
                
        % let's store the field name and the data to the responsePatterns
        % variable now
if all(all(isnan(beta_VC))) % we only store the data when there is data
    continue;
end

responsePatterns.(curr_epoch)(kk).name = fn; 
responsePatterns.(curr_epoch)(kk).rawdata = beta_VC; 
%responsePatterns.(curr_epoch)(kk).data = beta_VRC_has_data; % data without zeros
kk = kk+1;

    end
end
end
end

% %%
% for ii = 1:length(responsePatterns_delay)
% if isempty(responsePatterns_delay(ii).data)
%     responsePatterns_delay(ii) = [];
%     responsePatterns_stim(ii) = [];
% end
% end

%%
save(fullfile(out_dir,'responsePattern.mat'),'responsePatterns','good_cond')
end


%% Prepare model RDMs
% Prepare the model RDMs.
% clear
% clc
out_dir = rootDir;%'/mnt/Data1/linjdata1/vswmda/scan_data/rsa/beta_run_rsa_censor/'; % output directory
load(fullfile(out_dir,'f17','responsePattern.mat'),'good_cond') % only load the condition labels
good_cond([5:8,13:16]) = [];


good_cond = cellfun(@(C) C(1:3),good_cond,'UniformOutput',false)

% Recall eight patterns, corresponding to 8 columns of the 'cond' variable
% as well as the response matrix constructed above
% 3WL: 3dva, wide, L->R (0,0,0)
% 3TL: 3dva, tall, L->R (0,1,0)
% 3WR: 3dva, wide, R->L (0,0,1)
% 3TR: 3dva, tall, R->L (0,1,1)
% 5.5WL: 5.5dva, wide, L->R (1,0,0)
% 5.5TL: 5.5dva, tall, L->R (1,1,0)
% 5.5WR: 5.5dva, wide, R->L (1,0,1)
% 5.5TR: 5.5dva, tall, R->L (1,1,1)

% we will construct four different models:
% first model: UniqPat
% All eight patterns are unique
model.name{1} = 'UniqPat';
model.data{1} = double(~eye(8));
% It's important to remember that we will compare the model to the
% dissimilarity matrix. 0 is completely similar while 1 is completely
% dissimilar. So for the unique pattern model, we want diagonal components
% to be all zeros and all other components to be ones.

% second model: spatiotemporal pattern
% only (spatio)temporal patterns differ, such that left-to-right differs
% from right-to-left
model.name{2} = 'TempPat';
model.data{2} = repmat([zeros(2,2) ones(2,2);ones(2,2) zeros(2,2)],2,2);

% third model: spatial (shape) pattern
% only spatial-shape patterns differ, such that wide patterns differs
% from tall patterns
model.name{3} = 'ShapePat';
model.data{3} = double(repmat(~eye(2),4,4));

% fourth model: spatial set pattern
% only spatial-eccentricity patterns differ, such that 3-dva patterns differs
% from 5.5-dva patterns
% model.name{4} = 'SetPat';
% model.data{4} = repmat([zeros(4,4) ones(4,4);...
%                  ones(4,4) zeros(4,4)],4,4);
model.name{4} = 'RegPat';
model.data{4} = [zeros(4,4) ones(4,4);...
                  ones(4,4) zeros(4,4)];
             
% we want to plot those models and see if it makes sense
figure(1);clf
for pp = 1:4
subplot(2,2,pp)
heatmap(model.data{pp})
ax = gca;
ax.XData = good_cond;
ax.YData = good_cond;
ax.title(model.name{pp})
end

save(fullfile(out_dir,'model_exp2.mat'),'model','good_cond')
saveas(figure(1),fullfile(out_dir,'model_exp2.jpg'))

end