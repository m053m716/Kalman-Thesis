function fig = cometKinematicTrial(P,marker,iTrial)
%COMETFILTER Make comet tail for kinematic marker
%
%  fig = preview.cometFilter(wp,marker,iTrial);
%
% Inputs
%  P - Table of kinematic data, after cleaning step
%
% Output
%  fig - Figure handle
%
% See also: cpl

if ~strcmp(P.Properties.UserData.status,'clean')
   error('Need to do cleaning step first (see kal.getCleanKinematics)');
end

if nargin < 2
   marker = 'd1_d';
end

if nargin < 3
   iTrial = randi(max(P.TrialID),1);
end

xp = P.Xpred(P.TrialID==iTrial & strcmpi(P.Marker,marker) & strcmpi(P.Dimension,'X'),:);
yp = P.Xpred(P.TrialID==iTrial & strcmpi(P.Marker,marker) & strcmpi(P.Dimension,'Y'),:);
zp = P.Xpred(P.TrialID==iTrial & strcmpi(P.Marker,marker) & strcmpi(P.Dimension,'Z'),:);

xo = P.X(P.TrialID==iTrial & strcmpi(P.Marker,marker) & strcmpi(P.Dimension,'X'),:);
yo = P.X(P.TrialID==iTrial & strcmpi(P.Marker,marker) & strcmpi(P.Dimension,'Y'),:);
zo = P.X(P.TrialID==iTrial & strcmpi(P.Marker,marker) & strcmpi(P.Dimension,'Z'),:);

t = P.Properties.UserData.t;
tq = linspace(min(t),max(t),1000);

xp = interp1(t,xp,tq);
yp = interp1(t,yp,tq);
zp = interp1(t,zp,tq);

xo = interp1(t,xo,tq);
yo = interp1(t,yo,tq);
zo = interp1(t,zo,tq);

fig = figure('Name',sprintf('Kinematic Trial: %s - %d',marker,iTrial(1)),...
   'Color','w','Units','Normalized','Position',[0.2 0.2 0.4 0.4]); 
ax = axes(fig,'XColor','k','YColor','k','ZColor','k','View',[-37.5 30],...
   'FontName','Arial','NextPlot','add',...
   'XLim',[min([xo,xp]) max([xo,xp])],'YLim',[min([yo,yp]) max([yo,yp])],'ZLim',[min([zo,zp]) max([zo,zp])]);
box(ax,'on');
grid(ax,'on');
xlabel(ax,'X (a.u.)','FontName','Arial','Color','k');
ylabel(ax,'Y (a.u.)','FontName','Arial','Color','k');
zlabel(ax,'Z (a.u.)','FontName','Arial','Color','k');
title(ax,['Observed' newline getLabel(marker,iTrial)],...
   'FontName','Arial','Color','k');
comet3(ax,xo,yo,zo);
title(ax,['Predicted' newline getLabel(marker,iTrial)],...
   'FontName','Arial','Color','k');
comet3(ax,xp,yp,zp);
function str = getLabel(marker,iTrial)
   str = sprintf('Kinematic Trial: %s - %d',strrep(marker,'_','\_'),iTrial);
end

end