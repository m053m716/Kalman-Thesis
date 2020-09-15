function fig = crossTrialLFP(L,probe,channel,varargin)
%CROSSTRIALLFP Preview median of LFP across trials of same channel
%
%  fig = preview.crossTrialLFP(L,probe,channel,'Name',value,...)
%
% Inputs
%  L - LFP table as returned by L = init.lfpData();
%
% Output
%  fig - Figure handle
%
% Show trial-averaged values for LFP on a given channel
%
% See also: preview, utils, utils.getKinematicStates, init, init.kinData

pars = struct;
pars.AutoSave = false;
pars.CrossTrialFigureName = 'Cross-Trial LFP %s Ch-%s';
pars.XLim = [-0.25 0.25];
pars.XLabel = 'Time (sec)';
pars.XTick = [-0.25 -0.05 0 0.05 0.25];
pars.YLim = [-5 5];
fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

if nargin < 2
   up = unique(L.Probe);
   probe = up(randi(numel(up),1));
end

if nargin < 3
   uch = unique(L.Channel);
   channel = uch(randi(numel(uch),1));
end

fig = figure('Name',sprintf('Preview of Cross-Trial LFP (%s Ch-%s)',probe,channel),...
   'Color','w','Units','Normalized','Position',[1.085 0.084 0.322 0.336]); 

ttxt = sprintf('%s Ch-%s',probe,channel);
ax = axes(fig,...
      'XColor','k','YColor','k',...
      'XGrid','on','XLim',pars.XLim,'XTick',pars.XTick,...
      'LineWidth',1.5,'FontName','Arial',...
      'YLim',pars.YLim,...
      'NextPlot','add');
title(ax,ttxt,'FontName','Arial','Color','k'); 
ylabel(ax,'Z-Score','FontName','Arial','Color','k'); 
xlabel(ax,pars.XLabel,'FontName','Arial','Color','k');

idx = strcmpi(L.Probe,probe) & strcmpi(L.Channel,channel);
X = nanmedian(L.LFP(idx,:),1);
CI = iqr(L.LFP(idx,:),1);

t = L.Properties.UserData.t;
xcb = [t, fliplr(t)];
ycb = [X+(CI./2), fliplr(X)-(fliplr(CI./2))];

patch(ax,xcb,ycb,[0.5 0.5 0.5],...
   'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5,'EdgeColor','none',...
   'DisplayName','IQR');
line(ax,t,X,...
   'Color',[0.35 0.35 0.35],'LineWidth',1.5,'DisplayName','median');
legend(ax,'FontName','Arial','TextColor','black','EdgeColor','none','Color','none');

if ~pars.AutoSave
   return;   
end

pause(0.5);

if exist('figures/LFP','dir')==0
   mkdir('figures/LFP');
end

str = sprintf(pars.CrossTrialFigureName,probe,channel);
savefig(fig,fullfile('figures','LFP',[str '.fig']));
saveas(fig,fullfile('figures','LFP',[str '.png']));
delete(fig);
disp('Figure saved.');

end