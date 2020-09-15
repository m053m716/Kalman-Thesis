function SPK = loadSortedSpikes(varargin)
%LOADSORTEDSPIKES  Load (sorted) spike trains from BLOCK structure
%
%  SPK = cpl.loadSortedSpikes;
%  SPK = cpl.loadSortedSpikes('NAME',value,...);
%
%  --------
%   INPUTS
%  --------
%  varargin       :     (Optional) 'NAME', value input argument pairs.
%
%  --------
%   OUTPUT
%  --------
%   SPK           :     Formatted table with columns corresponding to X and
%                          fs, as well as additional columns corresponding
%                          to spike filename (and metadata?)

%% DEFAULTS
X = nan;
FS = 24414.0625;
DIR = 'D:\MATLAB\Projects\Data\R18-68_2018_07_24_1';

DEF_DIR = 'D:\MATLAB\Projects\Data\R18-68_2018_07_24_1';
SPIKE_DIR = '_wav-sneo_CAR_Spikes';
CLUS_DIR = '_wav-sneo_SPC_CAR_Sorted';
SPIKE_ID = 'ptrain';
CLUS_ID = 'sort';

DEBUG = false;

%% PARSE VARARGIN
for iV = 1:2:numel(varargin)
   eval([upper(varargin{iV}) '=varargin{iV+1};']);
end

%% SELECT BLOCK IF NECESSARY
if isnan(DIR)
   DIR = uigetdir(DEF_DIR,'Select recording BLOCK');
   if DIR==0
      error('No block selected, script aborted.');
   end
else
   if exist(DIR,'dir')==0
      error('Invalid BLOCK location. Check DIR (%s)',DIR);
   end
end

block = strsplit(DIR,filesep);
block = block{end};

%% LOAD ALL SPIKE TRAINS INTO CELL ARRAY AND GET SAMPLING RATE IF PRESENT
F = dir(fullfile(DIR,[block SPIKE_DIR],['*' SPIKE_ID '*.mat']));
C = dir(fullfile(DIR,[block CLUS_DIR],['*' CLUS_ID '*.mat']));

if ~iscell(X)
   X = cell(numel(F),1);
   fprintf(1,'\nLoading spike trains for %s...',block)
   h = waitbar(0,'Please wait, loading spike trains...');
   for ii = 1:numel(F)
      in = load(fullfile(F(ii).folder,F(ii).name));  
      cl = load(fullfile(C(ii).folder,C(ii).name));
      if DEBUG
         mtb(X); %#ok<UNRCH>
      end
      if ii == 1
         if isfield(in,'pars')
            if isfield(in.pars,'FS')
               fs = in.pars.FS;
               fprintf(1,'\n->\tFS detected: %g Hz\t\t\t\t ...',FS);
            else
               beep;
               fprintf(1,'\n->\tUsing default FS: %g Hz\t\t\t\t ...',FS);
               fs = FS;
            end
         else
            beep;
            fprintf(1,'\n->\tUsing default FS: %g Hz\t\t\t\t ...',FS);
            fs = FS;
         end
      end
      X{ii} = find(in.peak_train)./fs;
      X{ii}(cl.class==1) = [];      
      waitbar(ii/numel(F));
   end
   delete(h);
%    if exist('alert.mat','file')==0
%       beep;
%    else
%       alertSound = load('alert.mat','fs','sfx');
%       s = alertSound.sfx./max(abs(alertSound.sfx));
%       sound(s.*db2mag(-20),alertSound.fs);
%    end
   fprintf(1,'complete.\n');
else
%    if exist('alert.mat','file')==0
%       beep;
%       pause(0.5);
%       beep;
%    else
%       alertSound = load('alert.mat','fs','sfx');
%       s = alertSound.sfx./max(abs(alertSound.sfx));
%       sound(s.*db2mag(-20),alertSound.fs);
%       pause(0.5);
%       sound(s.*db2mag(-20),alertSound.fs);
%    end
   fprintf(1,'\n->\tUsing default FS: %g Hz\t\t\t\t ...\n',FS);
   fs = FS;
end
%% FORMAT TO TABLE
fname = {F.name}.';

ch = string(fname);
ch = strrep(extractBetween(ch,28,36),'_','');
probe = extractBefore(ch,'Ch');
chan = extractAfter(ch,'Ch');
SPK = table(probe,chan,X,'VariableNames',{'Probe','Channel','Spikes'});
SPK.Properties.UserData = struct('fs',fs);

end