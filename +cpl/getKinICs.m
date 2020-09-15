function [worldZ,Zica,W,T,mu] = getKinICs(worldPoints,varargin)
%% GETKINICS   Get Kinematic independent components

%% DEFAULTS
VAR_IDX = [1,3,4,6];
R = 3;
DO_ICA = true;

%% PARSE VARARGIN
for iV = 1:2:numel(varargin)
   eval([upper(varargin{iV}) '=varargin{iV+1};']);
end

%%
nTrials = numel(worldPoints);
nSamples = size(worldPoints{1}{1},1);

Z = nan(numel(VAR_IDX)*3,nSamples*nTrials);

%% Organize data
for ii = 1:nTrials
   vec = (nSamples*(ii-1)+1):(nSamples*ii);
   varIdx = 1:3;
   for iVar = 1:numel(VAR_IDX)
      v = VAR_IDX(iVar);
      Z(varIdx,vec) = worldPoints{ii}{v}.';      
      varIdx = varIdx + 3;
   end
end

%% Do kICA
[Zica,W,T,mu] = stat.kICA(Z,R);

%% Put data into similar structure
worldZ = cell(size(worldPoints));
for ii = 1:nTrials
   vec = (nSamples*(ii-1)+1):(nSamples*ii);
   if DO_ICA
      worldZ{ii} = Zica(:,vec).';
   else
      worldZ{ii} = Z(:,vec).';
   end
end

end