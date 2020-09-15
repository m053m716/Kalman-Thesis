function Z = neuralData(S,L)
%NEURALDATA Initialize neural data that combines L (LFP table) and S (Spikes table)
%
%  Z = init.neuralData(S,L);
%
% Inputs
%  S - Spikes table
%  L - LFP table
%
% Output
%  Z - "posterior" (measurements/observations variable) -- "noisy" neural data
%
% See also: init, init.mainData, init.spikeData, init.lfpData

Z = outerjoin(S,L,...
   'Type','Left',...
   'Keys',{'TrialID','Probe','Channel'},...
   'MergeKeys',true,...
   'LeftVariables',setdiff(S.Properties.VariableNames,{'LFP','XC_sensory'}),...
   'RightVariables',{'LFP','XC_sensory'});
Z.Properties.Description = 'Neural Data Table';
Z = movevars(Z,'TrialID','before',1);
Z = movevars(Z,{'Trial','Reach'},'before','Grasp');
Z.Properties.UserData.type = 'Neural';
Z.Properties.VariableDescriptions{'LFP'} = L.Properties.VariableDescriptions{'LFP'};
Z.Properties.VariableUnits{'LFP'} = L.Properties.VariableUnits{'LFP'};
Z.Properties.VariableUnits{'XC_sensory'} = '\muV^2';
Z.Properties.VariableDescriptions{'XC_sensory'} = L.Properties.VariableUnits{'XC_sensory'};
Z.Properties.VariableDescriptions{'Spikes'} =S.Properties.VariableDescriptions{'Spikes'};
Z.Properties.VariableUnits{'Spikes'} = S.Properties.VariableUnits{'Spikes'};
Z.Properties.UserData.t_spikes = S.Properties.UserData.t;
Z.Properties.UserData.t_lfp = L.Properties.UserData.t;
Z.LB_fit(isnan(Z.LB_fit)) = 0;
Z.UB_fit(isnan(Z.UB_fit)) = 0;
Z.N_fit(isnan(Z.N_fit)) = 0;

[G,ID] = findgroups(Z.ChannelID);
iBad = splitapply(@(x)any(x),Z.Exclude,G);
for ii = 1:numel(ID)
   if iBad(ii)
      Z(Z.ChannelID==ID(ii),:) = [];
   end
end
end