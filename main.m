%MAIN Outline for repository workflow

% Clear workspace and load data if needed
clc;
clearvars -except P Z

if exist('P','var')==0
   P = getfield(load('KinematicData.mat','P'),'P');
end
if exist('Z','var')==0
   Z = getfield(load('NeuralData.mat','Z'),'Z');
end

%% Merge data to get main data table, T
T = init.mainData('P',P,'Z',Z); % Even if not supplied, this function should parse the main table

% About these data:
%
% T has one row for each trial.
%
% All data are cell arrays, with each cell corresponding to a trial and
% each row within a cell corresponding to a single timestep. All timesteps
% for these data are synchronized so each data cell should have the same
% number of rows.
%
% T has 3 data variables:
%  1. Z_spike : "Measurement" : Square-root transformed spike counts
%  2. Z_lfp   : "Measurement" : Z-scored local field potentials (LFP)
%  3. X       :    "State"    : 3D-reconstruction of digit/hand markers
%
% T.Properties.UserData.t contains the relative time of each sample to the
% Grasp (synchronization event of interest), in units of seconds.
%
% T.Properties.UserData.fs is the sample rate (Hz)

