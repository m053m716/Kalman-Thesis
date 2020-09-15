function [Z,Xc,Xt,P] = getKinematicStates(P,varargin)
%GETKINEMATICSTATES Convert initial kinematic table P into states for cleaning
%
%  [Z,Xc,Xt,Xd] = utils.getKinematicStates(P,'Name',value,...);
%
% Inputs
%  P - Table of kinematic data variables
%
% Output
%  [Z,Xc,Xt] - See kal.getCleanKinematics
%
% See also: kal, kal.getCleanKinematics

pars = struct;
pars.DataVariable = 'X';
pars.DoMedianSubtraction = true;
fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

Gz = findgroups(P.TrialID);
[Gd,CID] = findgroups(P(:,{'TrialID','Dimension'}));
Gc = findgroups(CID.Dimension);
Gt = findgroups(P(:,{'Marker','Dimension'}));
Xd = cell2mat(splitapply(@(X){nanmedian(X,1)},P.(pars.DataVariable),Gd));
Xc = cell2mat(splitapply(@(X){X(:)'},Xd,Gc))';
Xt = cell2mat(splitapply(@(X){nanmedian(X,1)},P.(pars.DataVariable),Gt))';
Z  = cell2mat(splitapply(@(Z){Z'},P.(pars.DataVariable),Gz));

if ~pars.DoMedianSubtraction
   return;
end

dimMedians = nanmedian(Xc,1);
markerMedians = nanmedian(Xt,1);
Xc = Xc - dimMedians;
Xt = Xt - markerMedians;
Z = Z - markerMedians;


if strcmpi(P.Properties.UserData.status,'dirty')
   % Replicate for each trial and subtract as well
   trialMedians = repmat(markerMedians',numel(unique(P.TrialID)),1);
   P.X = P.X - trialMedians;
   P.Properties.UserData.status = 'medians-removed';
end


end