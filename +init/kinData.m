function P = kinData(varargin)
%KINDATA Initialize kinematic data
%
%  P = init.kinData(varargin);

pars = struct;
pars.BehaviorFile = 'BehaviorData.mat';
pars.File = '3DData.mat';
pars.FS = 50; % Target (decimated) sample rate
pars.BinCenters = linspace(-0.3,0.3,91); % Parameters from how pars.WP was extracted
pars.Vars = {'d1_d','d1_m','d1_p','d2_d','d2_m','d2_p'};
pars.VarEquivs = {'d1_d-x';'d1_d-y';'d1_d-z';...
                  'd1_m-x';'d1_m-y';'d1_m-z';...
                  'd1_p-x';'d1_p-y';'d1_p-z';...
                  'd2_d-x';'d2_d-y';'d2_d-z';...
                  'd2_m-x';'d2_m-y';'d2_m-z';...
                  'd2_p-x';'d2_p-y';'d2_p-z'};
pars.WP = [];

fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

if isempty(pars.WP)
   pars.WP = getfield(load(pars.File,'pars.WP'),'pars.WP');
end

nTrials = size(pars.WP,1);
nVars = numel(pars.VarEquivs);

data = cell(nTrials,1);
rmVars = setdiff(fieldnames(pars.WP{1}),pars.Vars);

h = waitbar(0,'Converting format of wp cell array struct...');
% Convert first so that we can decimate to correct sample rate
for iTrial = 1:nTrials
   pars.WP{iTrial} = rmfield(pars.WP{iTrial},rmVars);
   tmp = struct2table(pars.WP{iTrial});
   data{iTrial} = table2array(tmp(:,:))';
   
   waitbar(iTrial/nTrials);
end
delete(h);

TrialID = nan(nTrials*nVars,1);
Outcome = false(nTrials*nVars,1);
Marker = strings(nTrials*nVars,1);
behaviorData = getfield(load(pars.BehaviorFile,'behaviorData'),'behaviorData');
[data,t] = cpl.doDecimation(data,pars.BinCenters,pars.FS);
nSamples = numel(t);
X = nan(nTrials*nVars,nSamples);

h = waitbar(0,'Rearranging 3D motion variables table...');
iRow = 0;
for iTrial = 1:nTrials
   for ii = 1:nVars
      iRow = iRow + 1;
      TrialID(iRow) = iTrial;
      Marker(iRow) = string(pars.VarEquivs{ii});
      Outcome(iRow) = behaviorData.Outcome(iTrial);
      X(iRow,:) = data{iTrial}(ii,:);
   end
   waitbar(iTrial/nTrials);
end
delete(h);

P = table(TrialID,Outcome,Marker,X);
P.Dimension = extractAfter(P.Marker,"-");
P.Marker = extractBefore(P.Marker,"-");
P = movevars(P,'Dimension','before','X');

P.Properties.UserData = struct(...
   'type','Kinematic','t',t,'t0',pars.BinCenters,'fs',round(1/nanmean(diff(t))),...
   'status','dirty');
pars = rmfield(pars,'WP');
P.Properties.VariableUnits{'X'} = 'mm';
P.Properties.VariableDescriptions{'X'} = '3D "world" coordinates from stereo-calibrated cameras in X-Y-Z plane';
P.Properties.UserData.pars = pars;
P.Properties.Description = '3D-reconstruction coordinate table';
end