function frameTimes = parseBehaviorTimes(behaviorData,VideoOffset,varargin)
%% PARSEBEHAVIORTIMES  Get frame TIMES around behavior of interest
%
%  frameTimes = PARSEBEHAVIORTIMES(behaviorData,VideoOffset);
%  frameTimes = PARSEBEHAVIORTIMES(behaviorData,VideoOffset,'NAME',value,...);
%
%  --------
%   INPUTS
%  --------
%
%  behaviorData   :     Table with columns for event times, outcome,
%                       forelimb use, etc. Times are in seconds, relative
%                       to onset of NEURAL RECORDING.
%
%  VideoOffset    :     Scalar (seconds) of how many seconds prior to the
%                       video recording the neural recording was started.
%                       So for example if you started the Intan, then 30
%                       seconds later started the video, this number is 30.
%                       -> Can also be specified as a vector, for multiple
%                          camera offsets that correspond to the same
%                          neural recording.
%                          
%
%  varargin       :     (Optional) 'NAME', value input argument pairs.
%
%  --------
%   OUTPUT
%  --------
%
% By: Max Murphy  v1.0  10/03/2018  Original Version (R2017b)

%% DEFAULTS
EVENT = {'Grasp'};
E_PRE = -0.3;
E_POST = 0.3;

%% PARSE VARARGIN
for iV = 1:2:numel(varargin)
   eval([upper(varargin{iV}) '=varargin{iV+1};']);
end

%% BEGIN
N = size(behaviorData,1);
M = numel(EVENT);
C = numel(VideoOffset);

frameTimes = struct;
vec = [E_PRE,E_POST];

for iM = 1:M
   frameTimes.(EVENT{iM}) = cell(N,C);
   for iC = 1:C
      evt = behaviorData.(EVENT{iM}) - VideoOffset(iC);
      for iN = 1:N
         frameTimes.(EVENT{iM}){iN,iC} = vec + evt(iN);
      end
   end
end

end