function fig = plotMeanSlowWaves(snips,varargin)
%% PLOTMEANSLOWWAVES    Plot average of slow oscillations (<= 5 Hz)
%
%  fig = PLOTMEANSLOWWAVES(snips);
%  fig = PLOTMEANSLOWWAVES(snips,'NAME',value,...);
%
%  --------
%   INPUTS
%  --------
%    snips     :     Table from GETTRIALWAVEFORMS.
%
%  varargin    :     (Optional) 'NAME', value input argument pairs.
%
%
%  --------
%   OUTPUT
%  --------
%    fig       :     Figure handle to figure that has subplots
%                    corresponding to each channel, with each plot showing
%                    a snippet of the average slow waveform around grasps.
%
% By: Max Murphy  v1.0  11/21/2018  Original version (R2017b)

%% DEFAULTS
E_PRE = 500;  % ms
E_POST = 250; % ms
FS = 20;      % kHz
Y_LIM = [-3.5 3.5]; % z-score (?)

LINE_WIDTH = 2;
LINE_COLOR = 'k';

%% PARSE VARARGIN
for iV = 1:2:numel(varargin)
   eval([upper(varargin{iV}) '=varargin{iV+1};']);
end

%%
fig = figure('Name','Average Slow Oscillations',...
       'Units','Normalized',...
       'Color','w',...
       'Position',[0.2 0.2 0.6 0.7]);

for ii = 1:size(snips,1)
   subplot(5,10,ii); 
   plot((-E_PRE:(1/FS):E_POST),mean(snips.Slow{ii}),...
      'Color','k','LineWidth',2);
   xlim([-E_PRE E_POST]); 
   ylim(Y_LIM);
   
   addLabels(snips,ii);
   
end

   % Add correct label to subplot
   function addLabels(snips,ii)
      if ii >= 41
         xlabel('Time (ms)','Color','k','FontName','Arial');
      end
      if rem(ii-1,10)==0
         ylabel('Z-Score','Color','k','FontName','Arial');
      end
      
      title([snips.Probe{ii} '-' snips.Channel{ii}],...
         'FontName','Arial','Color','k');
      
   end

end