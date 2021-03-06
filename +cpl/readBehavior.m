function behaviorData = readBehavior(fname,tEvent)
%READBEHAVIOR  Read behavioral spreadsheet
%
%  behaviorData = cpl.readBehavior(fname)
%
%  --------
%   INPUTS
%  --------
%    fname     :     Filename of xlsx file with reach onset and grasp onset
%                    frames.
%
%    tEvent    :     Struct from CPL_OPTISYNC with information about video
%                    start and stop times relative to recording.
%
%  --------
%   OUTPUT
%  --------
%  behaviorData   :  Table with reach onset, grasp onset, outcome, and hand
%                    used.
%
% By: Max Murphy  v1.0  05/03/2018  Original version (R2017b)

% IMPORT DATA AND FORMAT
behaviorData = readtable(fname,'ReadVariableNames',true);
behaviorData(isnan(behaviorData.Reach),:) = [];

behaviorData.Button = [];

if ~isinf(tEvent.sync.stop)
   tVid = linspace(tEvent.sync.start,tEvent.sync.stop,tEvent.sync.nFrames);
else
   tVid = tEvent.sync.start:...
         (1/tEvent.sync.frameRate):...
         ((tEvent.sync.nFrames-1)/tEvent.sync.frameRate);
end
      
behaviorData.Reach = tVid(behaviorData.Reach).';
behaviorData.Grasp = tVid(behaviorData.Grasp).';

% CHECK IF FILE EXISTS, SAVE IF NOT
fname_out = strrep(fname,'.xlsx','.mat');
if exist(fname_out,'file')==0
   save(fname_out,'behaviorData','-v7.3');
end


end