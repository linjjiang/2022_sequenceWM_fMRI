% This script is to get beta map labels from the outputs of first-level GLM
% in SPM
% Beta map labels: which condition the map belongs to, and which run it
% belongs to
% We need these labels to extract specific beta maps for subsequent
% multivariate pattern analysis 

% By: Linjing Jiang
% Date: 03/26/2024

% Input:
% Which subject do you want to analyze: sb
% First-level GLM results from SPM (SPM.mat) in the corresponding folder

% Output:
% stim_beta: contains run_index - run number, condition_name - condition
% name, and beta_index - index of the beta map for the corresponding run
% and condition, for stimulus 
% delay_beta: same as stim_beta, except for delay period
% resp_beta: same as stim_beta, except for response period

% %% clean up the environment
% clear 
% close all
% clc

function get_beta_labels(rootDir,idx)
%% which subject, and which session do you want to analyze?

% data directory
%rootDir = '/mnt/Data1/linjdata1/vswmda/scan_data/spm/output/preproc_st_spatiotemporal_censor'; % root directory is the current directory
subDir = dir(fullfile(rootDir,'f*')); % subject directory. change if necessary.

% task
task_id = 'smgs'; % task name

for sb = idx%[1,2,3,6,7]%1:length(subDir) % for each subject %
    
    sub_id = subDir(sb).name; % subject id
    sessDir = dir(fullfile(subDir(sb).folder,subDir(sb).name,'S*'));
    
    for ss = 1:length(sessDir) % for each session under that subject
        sess_id = sessDir(ss).name; 
        
        clearvars -except rootDir subDir sessDir sb ss sub_id sess_id task_id idx
       
             % output directory
outputDir = fullfile(rootDir,sub_id,sess_id,task_id); % spm mat output directory
spmfile = fullfile(outputDir,'SPM.mat');
if ~isempty(dir(spmfile)) % if the output directory is not empty (there is already a SPM result)
    % we load the spm results file
    load(spmfile)
    
    % get the condition label
    for cc = 1:size(SPM.Vbeta,2)
    cond_label{cc} = SPM.Vbeta(cc).descrip;
    end
    
    % string of 8 conditions
    cond = {'_3WL','_3TL','_3WR','_3TR','_5.5WL','_5.5TL','_5.5WR','_5.5TR'};
    
    % find the beta map index of all 8 conditions, stimulus period
    % actual condition names
    cond_names = cellfun(@(C) ['stim' C],cond,'UniformOutput',false);
    stim_cond_idx = [];
    stim_cond_names = {};
    stim_run_idx = [];
    stim_beta_names = {};
    clc % clear command window
    for cc = 1:length(cond)  
        clear beta_name cond_name run_idx run_name temp
        
        temp = find(contains(cond_label,cond_names{cc}))';
        if ~isempty(temp)
    stim_cond_idx = [stim_cond_idx;temp]; % find beta map index of each condition
    stim_cond_names = [stim_cond_names;repmat(cond_names(cc),length(temp),1)];
    
    for tt = 1:length(temp)
        % get run number
    cond_name = split(cond_label{temp(tt)},' '); run_name = cond_name{4}; 
    run_idx(tt,1) = str2num(erase(run_name,{'Sn','(',')'}));
    
    % get beta map number
    
    beta_name{tt,1} = ['beta_' erase(cond_name{2},{'(',')'}) '.nii'];
    
    fprintf('%s\n',cond_label{temp(tt)})
    end

    stim_run_idx = [stim_run_idx;run_idx];
    stim_beta_names = [stim_beta_names;beta_name];
    
        end
    end
    
    % we want to create a table and attach the corresponding names for each
    % row and column
    stim_beta = table(stim_run_idx,stim_cond_names,stim_cond_idx,stim_beta_names,'VariableNames',{'run_index','condition_name','beta_index','beta_name'});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the beta map index of all 8 conditions, delay period
    cond_names = cellfun(@(C) ['delay' C],cond,'UniformOutput',false);
    delay_cond_idx = [];
    delay_cond_names = {};
    delay_run_idx = [];
    delay_beta_names = {};
    clc % clear command window
    for cc = 1:length(cond)  
        clear beta_name cond_name run_idx run_name temp
        
        temp = find(contains(cond_label,cond_names{cc}))';
        if ~isempty(temp)
    delay_cond_idx = [delay_cond_idx;temp]; % find beta map index of each condition
    delay_cond_names = [delay_cond_names;repmat(cond_names(cc),length(temp),1)];
    
    for tt = 1:length(temp)
        % get run number
    cond_name = split(cond_label{temp(tt)},' '); run_name = cond_name{4}; 
    run_idx(tt,1) = str2num(erase(run_name,{'Sn','(',')'}));
    
    % get beta map number
    
    beta_name{tt,1} = ['beta_' erase(cond_name{2},{'(',')'}) '.nii'];
    
    fprintf('%s\n',cond_label{temp(tt)})
    end

    delay_run_idx = [delay_run_idx;run_idx];
    delay_beta_names = [delay_beta_names;beta_name];
   
        end
    end
    
    % we want to create a table and attach the corresponding names for each
    % row and column
    delay_beta = table(delay_run_idx,delay_cond_names,delay_cond_idx,delay_beta_names,'VariableNames',{'run_index','condition_name','beta_index','beta_name'});
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find the beta map index of all 8 conditions, response period
    cond_names = cellfun(@(C) ['response' C],cond,'UniformOutput',false);
    response_cond_idx = [];
    response_cond_names = {};
    response_run_idx = [];
    response_beta_names = {};
    clc % clear command window
    for cc = 1:length(cond)  
        clear beta_name cond_name run_idx run_name temp
        
        temp = find(contains(cond_label,cond_names{cc}))';
        if ~isempty(temp)
    response_cond_idx = [response_cond_idx;temp]; % find beta map index of each condition
    response_cond_names = [response_cond_names;repmat(cond_names(cc),length(temp),1)];
    
    for tt = 1:length(temp)
        % get run number
    cond_name = split(cond_label{temp(tt)},' '); run_name = cond_name{4}; 
    run_idx(tt,1) = str2num(erase(run_name,{'Sn','(',')'}));
    
    % get beta map number
    
    beta_name{tt,1} = ['beta_' erase(cond_name{2},{'(',')'}) '.nii'];
    
        fprintf('%s\n',cond_label{temp(tt)})
    end

    response_run_idx = [response_run_idx;run_idx];
    response_beta_names = [response_beta_names;beta_name];
    
        end
    end
    
    % we want to create a table and attach the corresponding names for each
    % row and column
    response_beta = table(response_run_idx,response_cond_names,response_cond_idx,response_beta_names,'VariableNames',{'run_index','condition_name','beta_index','beta_name'});
  
% save for each session
save(fullfile(outputDir,'beta_label.mat'),'stim_beta','delay_beta','response_beta')
end 
    end
end
end