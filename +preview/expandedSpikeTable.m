function fig = expandedSpikeTable(X,probe,channel,trialID,varargin)
%EXPANDEDSPIKETABLE Preview trial(s) spike rate estimates from statistical model in "expanded" rates table
%
%  fig = preview.expandedSpikeTable(X,probe,channel,trialID,'Name',value,...);
%  
% Inputs
%  X - "Expanded" rates table (see: stat.getCleanSpikeRates)
%  probe - Name of probe
%  channel - Name of channel
%  trialID - Trial or trials to plot
%
% Output
%  fig - Figure handle
%
% See also: preview, stat, stat.getCleanSpikeRates, init, init.spikeData


pars = struct;
pars.Knots = 1:6:60;
pars.PosteriorMarkerSize = 8;
pars.PosteriorMarkerColor = [0.2 0.2 0.8];
pars.PosteriorMarkerAlpha = 0.8;

pars.PriorMarkerSize = 8;
pars.PriorMarkerColor = [0.5 0.5 0.5];
pars.PriorMarkerAlpha = 0.7;

pars.SGOrder = 3;
pars.SGFrameLen = 21;
fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

if nargin < 4
   trialID = randi(max(X.TrialID),1);
end

plotName = sprintf('%s - %s De-Noised Spike Rates',probe,channel);

fig = figure(...
   'Name',plotName,'Color','w','Units','Normalized','Position',[0.2 0.2 0.5 0.5]);
ax = axes(fig,'XColor','k','YColor','k','NextPlot','add',...
   'FontName','Arial','LineWidth',1.5,...
   'XLim',[min(X.t),max(X.t)]);


iAll = (strcmpi(X.Probe,probe)) & (strcmpi(X.Channel,channel));
x = X(iAll,:);
[g,tAll] = findgroups(x.t);
[tAll,iSort] = sort(tAll,'ascend');
tAll = reshape(tAll,1,numel(tAll));
ymu = splitapply(@(n)nanmean(n),x.Spikes,g);
ymu = reshape(ymu(iSort),1,numel(ymu));

ymu_to_knot = sgolayfilt(ymu,pars.SGOrder,pars.SGFrameLen);
ymu_to_knot = [max(round(ymu_to_knot(1)),0),max(ymu_to_knot(2:(end-1)),0),max(round(ymu_to_knot(end)),0)];
pars.Knots = union(pars.Knots,[1,numel(ymu_to_knot)]);

pp = csape(tAll(pars.Knots),ymu_to_knot(pars.Knots),'clamped');
sp = fnval(pp,tAll);

for ii = 1:numel(trialID)
   
  
   idx = (X.TrialID==trialID) & iAll;
   t = X.t(idx);

   errorbar(ax,t,X.N_fit(idx),X.CI_fit(idx,2)-X.N_fit(idx),X.N_fit(idx)-X.CI_fit(idx,1),...
      'LineStyle','none','LineWidth',1.25,'Color','r','DisplayName','95% CI: Prior');
   line(ax,t,X.N_Spline(idx),'LineWidth',1.5,'Color','k','Marker','*',...
      'DisplayName','Spline Fit');
   line(ax,t,X.Spikes(idx),'LineWidth',1.5,'Color',[0.5 0.5 0.5],'Marker','d','DisplayName','Observed');
   
   
   line(ax,tAll,ymu,'LineWidth',2,...
      'Color','b','DisplayName','Trial-Average');
   line(ax,tAll,sp,'LineWidth',2,...
      'Color','m','LineStyle',':','DisplayName','Clamped Spline',...
      'Marker','s','MarkerIndices',pars.Knots,'MarkerFaceColor','b');
   
   
end
legend(ax,'TextColor','black','FontName','Arial','Location','best');
ylabel(ax,'Spikes/sec','FontName','Arial','Color','k');

end