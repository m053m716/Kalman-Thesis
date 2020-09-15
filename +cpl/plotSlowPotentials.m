function D = plotSlowPotentials(behaviorData,varargin)
%PLOTSLOWPOTENTIALS  Plot average RAW (or filtered) LFP.
%
%  cpl.plotSlowPotentials(behaviorData);
%  D = cpl.plotSlowPotentials(behaviorData,'NAME',value,...);
%
%  --------
%   INPUTS
%  --------
%  behaviorData   :  From cpl.readBehavior. Data table containing times and
%                                           trial outcomes / grouping
%                                           variable for trials.
%
%  --------
%   OUTPUT
%  --------
%  D - [2,2] cell array, top-left is reach left, bottom right is
%              grasp right
%
%  Plots out ensemble-averaged slow potentials for each channel and
%  alignment condition.
%
% See also: cpl, cpl.readBehavior

% DEFAULTS
DIR = nan;

DEF_DIR = 'P:\Rat\BilateralReach\Murphy';
RAW_DIR = '_RawData';
RAW_ID = '_Raw_';

E_PRE = 1.0;
E_POST = 0.5;

YLIM = [-50 50];
XLIM = [-E_PRE E_POST];
LOWPASS = 10;  % Hz

% PARSE VARARGIN
for iV = 1:2:numel(varargin)
   eval([upper(varargin{iV}) '=varargin{iV+1};']);
end

% GET DIRECTORY
if isnan(DIR)
   DIR = uigetdir(DEF_DIR,'Select recording BLOCK');
   if DIR == 0
      error('No recording BLOCK specified. Script aborted.');
   end
end

block = strsplit(DIR,filesep);
block = block{end};
F = dir(fullfile(DIR,[block RAW_DIR],['*' RAW_ID '*.mat']));

% GET INDEXING FOR BEHAVIOR TYPES
iLeft = ismember(behaviorData.Forelimb,'L') | ismember(behaviorData.Forelimb,0);
iRight = ismember(behaviorData.Forelimb,'R')| ismember(behaviorData.Forelimb,1);    
    
for ii = 1:numel(F)
   str_info = strsplit(F(ii).name(1:end-4),'_');
   ch = str_info{end};
   probe = str_info{end-2};
   figure('Name',[block ': ' probe ' - Channel ' ch],...
          'Color','w',...
          'NumberTitle','off',...
          'WindowStyle','docked');
       
   in = load(fullfile(F(ii).folder,F(ii).name));
   [b,a] = butter(4,2*LOWPASS/in.fs);
   data = filtfilt(b,a,double(in.data));    
   subplot(2,2,1);
   [avg,t,err,D{1,1}{ii}] = ...
      cpl.getRawAverage(data,in.fs,behaviorData.Reach(iLeft),...
            'E_PRE',E_PRE,...
            'E_POST',E_POST);
   errorbar(t,avg,err);
   title('Reach - L');
   ylim(YLIM);
   xlim(XLIM);
   
   subplot(2,2,2);
   [avg,t,err,D{1,2}{ii}] = ...
      cpl.getRawAverage(data,in.fs,behaviorData.Reach(iRight),...
            'E_PRE',E_PRE,...
            'E_POST',E_POST);
   errorbar(t,avg,err);
   ylim(YLIM);
   xlim(XLIM);
   title('Reach - R');
   
   
   subplot(2,2,3);
   [avg,t,err,D{2,1}{ii}] = ...
      cpl.getRawAverage(data,in.fs,behaviorData.Grasp(iLeft),...
            'E_PRE',E_PRE,...
            'E_POST',E_POST);
   errorbar(t,avg,err);
   ylim(YLIM);
   xlim(XLIM);
   title('Grasp - L');
   
   subplot(2,2,4);
   [avg,t,err,D{2,2}{ii}] = ...
      cpl.getRawAverage(data,in.fs,behaviorData.Grasp(iRight),...
            'E_PRE',E_PRE,...
            'E_POST',E_POST);
   errorbar(t,avg,err);
   ylim(YLIM);
   xlim(XLIM);
   title('Grasp - R');
   
   suptitle([strrep(block,'_','-') ': ' probe ' - Channel ' ch]);
end

end