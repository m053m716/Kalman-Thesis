function fig = PETH(t,X)
%PETH Just make a bar graph in a figure
%
%  fig = preview.PETH(S);
%  fig = preview.PETH(t,X);
%
% Inputs
%  S - Table of spikes as returned by S = init.spikeData
%  -- or --
%
%  t - Bin centers for spike counts
%  X - Spike counts
%
% Output
%  fig - Figure handle
%
% See also: init, init.spikeData, preview

POS = [0.4 0.4 0.2 0.2];
J = [0.1 0.1 0.0 0.0];
pos = POS + randn(1,4).*J;

if istable(t)
   X = t.Spikes;
   t = t.Properties.UserData.t;
end

fig = figure(...
   'Name','Spike Count PETH',...
   'Color','w',...
   'Units','Normalized',...
   'Position',pos);
ax = axes(fig,'NextPlot','add',...
   'XColor','k','YColor','k',...
   'LineWidth',1.5,'FontName','Arial');
xlabel(ax,'Time (ms)','FontName','Arial','Color','k');
ylabel(ax,'Spike Counts','FontName','Arial','Color','k');
bar(ax,round(t.*1e3),nansum(X,1),1,'EdgeColor','none','FaceColor',rand(1,3));

end