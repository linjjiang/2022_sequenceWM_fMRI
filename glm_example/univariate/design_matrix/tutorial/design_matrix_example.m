% Design matrix for fMRI
% By: Linjing
% Update: 03/29/2023

% This is an example script to create design matrices for your experiment.
% Note: this script only creates experimental regressors. You need other
% scripts to create motion regressors.

%% Clean up the coding environment
clear % clear the workspace
close all % close all figures
clc % clean up the command window before you get started

%% Load the data
% Let's use my data as an example.
% In my experiment, we recorded subjects' eye movements when they were
% performing a memory-guided saccade task. Subjects saw an image briefly on
% the screen and moved their eyes to the remembered image location after a
% lengthened delay (8-10 sec). 
% Here let's load data from one run of mgs:
load('data_mgs1.mat')
clearvars -except edf

% 'edf' contains eye data and experimental parameters from one mgs run.

%% Make design matrices

% You need to construct three variables, names, durations and onsets, for
% your spm design matrix.

% It is always good to take a look at the overall BOLD signal patterns
% regardless of which types of trial it is. 

% In my case, I want to ask how the brain is activated during memory
% retention regardless of trial type. To address this, I construct one 
% experimental regressor, named 'delay', that indicates when the delay
% period of the task starts. 

names = {'delay'}; % name of the regressor. Make sure you use the exact naming here, i.e., 'names' not 'name'.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The next few lines are to get the onsets and durations of 'delay'. 
%%% You should use your own scripts to get the timing of your task. %%%%%%
% Retrieve scan start time
scan_start = edf.default_events.Messages.time(contains(edf.default_events.Messages.info,'ScanTrig'));  % In eyelink time
% Get the timing of each event across trials. I have 16 trials and 7
% events. Delay start is event 4 (4th column) and delay end is event 5 (5th
% column).
time = (edf.events.msg.time - scan_start)/1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the onset and duration of the delay event. Here I model the onsets as
% the middle of the delay period, and duration as 0.
onsets = {(time(:,4)'+time(:,5)')/2};
durations = {zeros(size(time(:,4)'))};

% Now, save names, onsets and durations to a .mat file 
save('design_matrix_mgs1.mat','names','onsets','durations')

% You can apply the same type of codes to construct design matrices for
% other runs if you have multiple scanning runs. 

%% Make design matrices: advanced
% Let's say I have multiple types of trials. I want to go beyond examining 
% the overall activity to compare and contrasts BOLD activation across 
% trial types.

% To do that, you need to figure out which types of trials you want to
% model and construct separate regressors for each.

% An example: I want to contrast brain activity when image is memorized on
% the right half of the visual field compared to that on the left half of
% the visual field. Then I will construct two regressors, one indicating
% the delay timing of the right trials and the other showing timing of the
% left trials.

% First, figure out which trials have left vs. right presentation of
% stimuli. This infomation is stored in the 'edf.param.probe_quad'.
probe_quad = edf.param.probe_quad(1:16); % the quadrant of the stimulus for each trial (I hae 16 trials)
left_trial_index = find(probe_quad == 2 | probe_quad == 3);
right_trial_index = find(probe_quad == 1 | probe_quad == 4);

% Then, use the corresponding trial index to retrieve timing of each type
% of trials respectively.
names = {'delay_left','delay_right'}; % name of the regressor
onsets = {(time(left_trial_index,4)'+time(left_trial_index,5)')/2, ...
    (time(right_trial_index,4)'+time(right_trial_index,5)')/2}; % delay onsets of left and right trials
durations = {zeros(size(time(:,4)')), zeros(size(time(:,4)'))}; % durations

% Now, save names, onsets and durations to a .mat file 
save('design_matrix_advanced_mgs1.mat','names','onsets','durations')

% Again, you can apply the same code to multiple runs!
