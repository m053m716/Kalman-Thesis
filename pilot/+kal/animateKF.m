function CData = animateKF(F,trialIndex,iter,times)
%ANIMATEKF Return animation stack for a given trial to export filter weights video output
%
%  CData = kal.animateKF(F,trialIndex,times);
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

if nargin < 2
   trialIndex = randi(max(F.Trial),1);
end

if nargin < 3
   iter = randi(max(F.Iteration),1);
end

if nargin < 4
   times = F.ts;
end
idx = F.Trial==trialIndex & ismember(F.ts,times) & ismember(F.Iteration,iter);
F = F(idx,:);
CData = cat(3,F.K{:});
if nargout > 0
   return;
end

fig = figure('Name','Kalman Animation','Color','w','Units','Normalized',...
   'Position',[0.3 0.3 0.3 0.3]);
ax = axes(fig,'XColor','none','YColor','none','Color','none',...
   'NextPlot','add');
lab = title(ax,sprintf('Trial %d (%5.2f ms)',trialIndex,F.ts(1)),...
   'FontName','Arial','Color','k');

C = softmax(CData(:,:,1)')';

im = imagesc(ax,C);
colorbar(ax);
for ii = 2:size(CData,3)
   pause(0.25);
   C = softmax(CData(:,:,ii)')';
   set(im,'CData',C);
   set(lab,'String',sprintf('Trial %d (%5.2f ms)',trialIndex,F.ts(ii)));
   drawnow;
end
delete(fig);

end