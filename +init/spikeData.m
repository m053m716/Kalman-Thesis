function S = spikeData(varargin)
%SPIKEDATA Initialize spike data
%
%  S = init.spikeData('Name',value,...);
%
% Inputs
%  varargin - Optional name,value pairs (see pars, below)
%
% Output
%  S - Table of trial-aligned per-channel, per-trial spiking
%        -> Alignment is to the Grasp.
%
% See also: init, init.kinData, init.lfpData

pars = struct;
pars.Alignment = 'Grasp';
pars.BehaviorFile = 'BehaviorData.mat';
pars.Bins = -0.610:0.020:0.610; % Bin edges (seconds)
pars.SpikeTimes = [];

fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

behaviorData = getfield(load(pars.BehaviorFile,'behaviorData'),'behaviorData');
nTrials = size(behaviorData,1);
% behaviorData.TrialID = (1:nTrials)';

if isempty(pars.SpikeTimes)
   pars.SpikeTimes = cpl.loadSortedSpikes();
end
nChannels = size(pars.SpikeTimes,1);
nBins = numel(pars.Bins)-1;

S = repmat(behaviorData,nChannels,1);
S.Probe = strings(nTrials*nChannels,1);
S.Channel = strings(nTrials*nChannels,1);
S.Spikes = zeros(nTrials*nChannels,nBins);

% iThis = 0:nChannels:(nChannels*(nTrials-1));
h = waitbar(0,'Please wait, getting individual trial spike counts...');
idx = 0;
for iTrial = 1:nTrials
%    iThis = iThis+1;
   for iChannel = 1:nChannels
%       idx = iThis(iTrial);
      idx = idx + 1;
      S.Probe(idx) = pars.SpikeTimes.Probe(iChannel);
      S.Channel(idx) = pars.SpikeTimes.Channel(iChannel);
      S.Spikes(idx,:) = histcounts(pars.SpikeTimes.Spikes{iChannel}-behaviorData.(pars.Alignment)(iTrial),pars.Bins);
      S.TrialID(idx) = iTrial;
   end
   waitbar(iChannel/nChannels);
end
delete(h);
S.Properties.VariableUnits{'Trial'} = 'sec';
S.Properties.VariableUnits{'Reach'} = 'sec';
S.Properties.VariableUnits{'Grasp'} = 'sec';
S.Properties.VariableUnits{'Spikes'} = 'spikes';
S.Properties.VariableDescriptions{'Spikes'} = 'Spike counts, centered around the Grasp event time';
S.Properties.UserData = struct('fs_orig',pars.SpikeTimes.Properties.UserData.fs,...
   't',pars.Bins(1:(end-1))+nanmean(diff(pars.Bins)./2),...
   'pars',pars,'type','Spikes');
S.Properties.UserData.fs = round(1/nanmean(diff(S.Properties.UserData.t)));

S = movevars(S,{'TrialID','Probe','Channel'},'Before',1);

end