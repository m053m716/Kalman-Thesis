function fig = cometKinematicTrial(wp,marker,iTrial)
%COMETKINEMATICTRIAL Make comet tail for kinematic marker
%
%  fig = cpl.cometKinematicTrial(wp,marker,iTrial);
%
% Inputs
%  wp - Loaded from reconstruction data in Data/Kalman-Thesis
%
% Output
%  fig - Figure handle
%
% See also: cpl

if nargin < 2
   marker = 'd1_d';
end

if nargin < 3
   iTrial = randi(numel(wp),1);
end

x = wp{iTrial(1)}.(marker)(:,1);
y = wp{iTrial(1)}.(marker)(:,2);
z = wp{iTrial(1)}.(marker)(:,3);
fig = figure('Name',sprintf('Kinematic Trial: %s - %d',marker,iTrial(1)),...
   'Color','w','Units','Normalized','Position',[0.2 0.2 0.4 0.4]); 
ax = axes(fig,'XColor','k','YColor','k','ZColor','k','View',[-37.5 30],...
   'FontName','Arial','NextPlot','add',...
   'XLim',[min(x) max(x)],'YLim',[min(y) max(y)],'ZLim',[min(z) max(z)]);
box(ax,'on');
grid(ax,'on');
xlabel(ax,'X (mm)','FontName','Arial','Color','k');
ylabel(ax,'Y (mm)','FontName','Arial','Color','k');
zlabel(ax,'Z (mm)','FontName','Arial','Color','k');
lab = title(ax,getLabel(marker,iTrial(1)),'FontName','Arial','Color','k');
for ii = 1:numel(iTrial)
   set(lab,'String',getLabel(marker,iTrial(ii)));
   x = wp{iTrial(ii)}.(marker)(:,1);
   y = wp{iTrial(ii)}.(marker)(:,2);
   z = wp{iTrial(ii)}.(marker)(:,3);
   comet3(ax,x,y,z);
end

   function str = getLabel(marker,iTrial)
      str = sprintf('Kinematic Trial: %s - %d',strrep(marker,'_','\_'),iTrial);
   end

end