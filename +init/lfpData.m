function L = lfpData(varargin)
%LFPDATA Initialize LFP data
%
%  L = init.lfpData('Name',value,...)
%
% Inputs
%  varargin - (Optional) 'Name',value input argument pairs
%
% Output
%  L
%
% See also: init

pars = struct;
pars.File = 'SnippetData.mat';
pars.FS = 50;
pars.LFPVar = 'SG'; % Savitzky-Golay low pass filtered field potentials
pars.XCProbe = "P1";
pars.XCChannels = ["024","025"];
pars.Snips = [];
pars.TCorr = 0.2; % Time for cross-correlation window (seconds)

fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

if isempty(pars.Snips)
   pars.Snips = getfield(load(pars.File,'snips'),'snips');
end

tic; fprintf(1,'Decimating lowpass-filtered LFP data...');
[data,t] =  cpl.doDecimation(pars.Snips.(pars.LFPVar),...
                             pars.Snips.Properties.UserData.t,...
                             pars.FS);
fprintf(1,'complete (%5.2f sec)\n',toc);

p = pars.Snips.Probe;
ch = pars.Snips.Channel;

% Get dimensions of Tables, data matrices
nChannels = numel(ch);
nTrials = size(data{1},1);
nSamples = numel(t);
nRows = nChannels * nTrials;

TrialID = nan(nRows,1);
Channel = strings(nRows,1);
Probe = strings(nRows,1);
LFP = nan(nRows,nSamples);
XC_sensory = nan(nRows,nSamples-(pars.TCorr*pars.FS)+1);
iRow = 0;

% Get channels to use mean for cross correlation
chLFP = find(ismember(p,pars.XCProbe) & ismember(ch,pars.XCChannels));
nChXC = numel(chLFP);
tic; fprintf(1,'Exporting LFP data to table...');
h = waitbar(0,'Exporting LFP data to table...');
for iTrial = 1:nTrials
   mu = zeros(nSamples,1);
   for iCh = 1:nChXC
      mu = mu + (data{iCh}(iTrial,:)' .* (1/nChXC));
   end
   for iCh = 1:nChannels
      iRow = iRow + 1;
      TrialID(iRow) = iTrial;
      Channel(iRow) = string(ch{iCh});
      Probe(iRow) = string(p{iCh});
      LFP(iRow,:) = data{iCh}(iTrial,:);
      [c0,tc] = utils.computeZeroLagCorrelation(mu,LFP(iRow,:)',...
         'ts',t,'T',pars.TCorr);
      XC_sensory(iRow,:) = c0';
   end
   waitbar(iTrial/nTrials);
end
fprintf(1,'complete (%5.2f sec)\n',toc);
delete(h);
L = table(TrialID,Probe,Channel,LFP,XC_sensory);
L.Properties.VariableUnits{'LFP'} = pars.Snips.Properties.VariableUnits{pars.LFPVar};
L.Properties.VariableDescriptions{'LFP'} = pars.Snips.Properties.VariableDescriptions{pars.LFPVar};
pars = rmfield(pars,'Snips');
L.Properties.UserData = struct('type','LFP','t',t,'tc',tc,'fs',pars.FS,'pars',pars);
L.Properties.VariableDescriptions{'XC_sensory'} = 'Cross-correlation with "sensory LFP" (average of selected channels, see pars)';
end