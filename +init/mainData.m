function [T,Z] = mainData(varargin)
%MAINDATA Initialize main data table
%
%  [T,Z] = init.mainData('Name',value,...);
%
% Inputs
%  varargin - (Optional) 'Name',value input argument pairs (see pars)
%
% Output
%  T - Main data table that synthesizes:
%        S = init.spikeData;
%        P = init.kinData;
%        L = init.LFPData;
%
% See also: init, init.spikeData, init.kinData, init.LFPData

pars = struct;
pars.FS = 50; % Determined empirically: decimated sample rate
pars.FS_INTERP = 150;
pars.L = [];
pars.LFile = 'LFPData.mat';
pars.P = [];
pars.PFile = 'KinematicData.mat';
pars.S = [];
pars.SFile = 'SpikeData.mat';
pars.Snips = [];
pars.SpikeTimes = [];
pars.TROI = [-0.240, 0.240];
pars.WP = [];
pars.Z = []; % Main neural data table

fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   else
      warning('%s is not a member of pars. Check inputs.',varargin{iV});
   end
end

if isempty(pars.Z)
   if isempty(pars.L)
      if exist(pars.LFile,'file')==0
         pars.L = init.lfpData('Snips',pars.Snips);
      else
         pars.L = getfield(load(pars.LFile,'L'),'L');
      end
   end

   if isempty(pars.S)
      if exist(pars.SFile,'file')==0
         pars.S = init.spikeData('SpikeTimes',pars.SpikeTimes);
      else
         pars.S = getfield(load(pars.SFile,'S'),'S');
      end
   end

   pars.Z = init.neuralData(pars.S,pars.L);
end
Z = pars.Z;
ts = Z.Properties.UserData.t_spikes;
iSpike = (ts > pars.TROI(1)) & (ts < pars.TROI(2));
tl = Z.Properties.UserData.t_lfp;
iLFP = (tl > pars.TROI(1)) & (tl < pars.TROI(2));

[G,T] = findgroups(pars.Z(:,'TrialID'));
T.Z_spike = splitapply(@(Z){Z(:,iSpike)'},Z.N_fit,G);
T.Z_lfp   = splitapply(@(Z){Z(:,iLFP)'},Z.LFP,G); 
T.Z_xc    = splitapply(@(Z){Z(:,iLFP)'},Z.XC_sensory,G);
T.Properties.VariableUnits{'Z_spike'} = pars.Z.Properties.VariableUnits{'Spikes'};
T.Properties.VariableUnits{'Z_lfp'} = pars.Z.Properties.VariableUnits{'LFP'};
T.Properties.VariableUnits{'Z_xc'} = '\muV^2';
% % % % Get kinematic data % % % %
if isempty(pars.P)
   if exist(pars.PFile,'file')==0
      pars.P = init.kinData('WP',pars.WP,...
                            'BinCenters',pars.Z.Properties.UserData.t);
   else
      pars.P = getfield(load(pars.PFile,'P'),'P');
   end
end

% Apply cleaning to P %
pars.P = kal.getCleanKinematics(pars.P,'MakeFigures',true,'TROI',pars.TROI);
to = pars.P.Properties.UserData.t; % Fixes and accounts for "TROI"
tq = linspace(to(1),to(end),numel(to)*round(pars.FS_INTERP/pars.FS));

T.Properties.UserData = struct('type','Main');
T.Properties.UserData.nTrials = size(T,1);
T.Properties.UserData.nSamples = numel(tq);
T.Properties.UserData.t = tq;
T.Properties.UserData.fs = pars.FS_INTERP;

% X: "hidden" (state variable) -- kinematic states to be estimated
G = findgroups(pars.P(:,'TrialID'));
T.X = splitapply(@(X){interp1(to,X',tq)},pars.P.Xpred,G);
T.Xo = splitapply(@(X){interp1(to,X',tq)},pars.P.X,G);
T.Z_spike = cellfun(@(Z)interp1(to,Z,tq),T.Z_spike,'UniformOutput',false);
T.Z_lfp = cellfun(@(Z)interp1(to,Z,tq),T.Z_lfp,'UniformOutput',false);
T.Z_xc = cellfun(@(Z)interp1(to,Z,tq),T.Z_xc,'UniformOutput',false);
T.Properties.VariableUnits{'X'} = 'mm';
T.Properties.VariableDescriptions{'X'} = 'Positions after Kalman filter cleaning';
T.Properties.VariableUnits{'Xo'} = 'mm';
T.Properties.VariableDescriptions{'Xo'} = 'Observed positions prior to cleaning';


[~,T.Properties.UserData.Channels] = findgroups(pars.Z(:,{'Probe','Channel'}));
T.Properties.UserData.Z = pars.Z;
T.Properties.UserData.P = pars.P;
pars = rmfield(pars,{'L','Snips','P','WP','S','SpikeTimes','Z'}); % Remove large data tables from parameters struct
T.Properties.UserData.pars = pars;
T.Properties.Description = 'Merged table of trials for Kalman Filter with variables split into State (X) and Measurement (Z) datasets';
end