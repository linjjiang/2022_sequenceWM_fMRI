% Extract 6 motion regressors from halfpipe outputs
% By: Linjing Jiang
% Update: 03/29/2023

% This script shows an example of how to extract motion regressors from
% halfpipe output. Before you get started, navigate to your halfpipe output
% folder 'derivatives/halfpipe/subject_id/func/'. The file we are going to 
% work on is called 'xxxxxxx-confounds_regressors.tsv', which contains a
% lot of different confounding variables output from preprocessing. We only
% need 6 of them: 
% trans_x: x transition
% trans_y: y transition
% trans_z: z transition
% rot_x: x rotation
% rot_y: y rotation
% rot_z: z rotation

% For tutorial demonstration, here I provided an example confounds
% regressor output from halfpipe, called 'mgs1_confounds_regressors.tsv'.
% You should find similar tsv files in your derivatives folder but with
% lengthenend names. 

% To get started, clean up the environment
clear
close all
clc

% Load the data
data = importdata('mgs1_confounds_regressors.tsv', '\t');

% Extract 6 motion parameters
% Find the column number corresponding to these 6 motion parameters
cols = ismember(data.textdata,{'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'});
% Extract these 6 columns
motion_regressor = data.data(:,cols);

% Save motion regressors to a .txt file
writematrix(motion_regressor, 'mgs1_motion_regressors.txt');

% You can apply the same code to different runs respectively!
