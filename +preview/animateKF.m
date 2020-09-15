function [fig,CData] = animateKF(K_k,trialIndex,varargin)
%ANIMATEKF Return animation stack for a given trial to export filter weights video output
%
%  fig = preview.animateKF(K_k,trialIndex,varargin);
%  [fig,CData] = preview.animateKF(K_k,trialIndex,varargin);
%
% Inputs
%  F  -  Table returned by kal.estimateKF
%  trialIndex - (Optional) trial index to export video for
%  times - Time range for video
%
% Output
%  CData - Image stack where third dimension is a new frame
%
% See also: kal, kal.estimateKF

pars = struct;
pars.CLim = [0 0.2];
pars.FrameRateApparent = 5;
% pars.XTick = [0.5, 17.5];
pars.XTick = 0.5:17.5;
% pars.XTickLabels = ["Kin^+_{1}","Kin^+_{-18}"];
pars.XTickLabels = ["Dig-1^X_d","Y","Z",...
                    "Dig-1^X_m","Y","Z",...
                    "Dig-1^X_p","Y","Z",...
                    "Dig-2^X_d","Y","Z",...
                    "Dig-2^X_m","Y","Z",...
                    "Dig-2^X_p","Y","Z"];
% pars.YTick = [0.5,24.5,25.5,26.5,28.5,29.5,31.5,48.5]; % States
% pars.YTickLabels = ["Spikes-01","Spikes-25","LFP-M1","LFP-S1","M1-S1-XC","S1-S1-XC","Kin^-_{1}","Kin^-_{-18}"];
pars.YTick = [0.5,24.5,25.5,26.5,28.5,29.5,30.5]; % States
pars.YTickLabels = ["Spikes-01","Spikes-25","LFP-M1","LFP-S1","M1-S1-XC","S1-S1-XC","S1-S1-XC"];
fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV+1});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

if nargin < 2
   trialIndex = randi(size(K_k,1),1);
end


if nargin < 3
   times = linspace(-0.22,0.22,69);
end

CData = K_k{trialIndex};
if nargout > 1
   return;
end

fig = figure(...
   'Name','Kalman Animation',...
   'Color','w','Units','Normalized',...
   'Position',[0.3 0.1 0.3 0.8]);
ax = axes(fig,...
   'XColor','k',...
   'TickDir','out',...
   'XLim',[0,pars.XTick(end)+0.5],...
   'XTick',pars.XTick,...
   'XTickLabels',pars.XTickLabels,...
   'XTickLabelRotation',60,...
   'YTick',pars.YTick,...
   'YTickLabels',pars.YTickLabels,...
   'YLim',[0,pars.YTick(end)+0.5],...
   'YColor','k',...
   'YDir','reverse',...
   'Color','none',...
   ...'CLim',pars.CLim,...
   'NextPlot','add');
xlabel(ax,'Measurement','FontName','Arial','Color','k');
ylabel(ax,'State','FontName','Arial','Color','k');
times = times.*1e3;
lab = title(ax,sprintf('Trial-%02d (%5.2f ms)',...
   trialIndex,times(1)),...
   'FontName','Arial','Color','k');

C = softmax(CData(:,:,1)')';

x = [pars.XTick(1),pars.XTick(end)];
y = [pars.YTick(1),pars.YTick(end)];
im = imagesc(ax,x,y,C);
% set(im,'CDataMapping','scaled');
colorbar(ax);
for ii = 2:size(CData,3)
   pause(1/pars.FrameRateApparent);
   C = softmax(CData(:,:,ii)')';
   set(im,'CData',C);
   set(lab,'String',sprintf('Trial-%02d (%+5.2f ms)',trialIndex,times(ii)));
   drawnow;
end
if nargout < 1
   delete(fig);
end

end