function byChannel = rasterByTrial2ByChannel(X)
%RASTERBYTRIAL2BYCHANNEL   Convert cell arrays from rasters by trial with all units to by unit for all trials
%
%  byChannel = cpl.rasterByTrial2ByChannel(X);
%
%  --------
%   INPUTS
%  --------
%     X     :     Cell array, where each cell contains binned rows that
%                 correspond to a single channel's activity, aligned to a
%                 single trial.
%
%  --------
%   OUTPUT
%  --------
%  byUnit   :     Cell array, where each cell contains binned rows that
%                 correspond to a single trial, and a given cell element
%                 corresponds to units of a single channel.

% PARSE INFO FROM INPUT CELL ARRAY
nTrial = numel(X);
nChannel = size(X{1},1);
nBin = size(X{1},2);

% LOOP THROUGH AND CHANGE CELL/ROW FORMAT
byChannel = cell(nChannel,1);
for iChannel = 1:nChannel
   byChannel{iChannel} = nan(nTrial,nBin);
   for iTrial = 1:nTrial
      byChannel{iChannel}(iTrial,:) = X{iTrial}(iChannel,:);
   end
end

end