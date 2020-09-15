function fig = medianKinematics(P,varargin)
%MEDIANKINEMATICS Preview median of kinematic components to do subtract
%
%  fig = preview.medianKinematics(P);
%  fig = preview.medianKinematics(P,'Name',value,...);
%
% Inputs
%  P - Kinematics table as returned by P = init.kinData();
%
% Output
%  fig - Figure handle
%
% Show the median values that will be used for "state estimator" cleaning
% using Kalman formulation for noisy data returned by AI.
%
% See also: preview, utils, utils.getKinematicStates, init, init.kinData

pars = struct;
pars.AutoSave = false;
pars.DataVariable = 'X';
pars.CrossTrialFigureName = 'Cross-Trial Medians for State Data';
pars.CrossStateFigureName = 'Cross-State Medians Across Trials';
pars.XLim = [-0.25 0.25];
pars.XLabel = 'Time (sec)';
pars.XTick = [-0.25 -0.05 0 0.05 0.25];
fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

if isa(P,'table')
   test = struct; 
   [test.Z,test.Xc,test.Xt] = utils.getKinematicStates(P,...
      'DoMedianSubtraction',false,...
      'DataVariable',pars.DataVariable);
   t = P.Properties.UserData.t;
else
   test = P;
   t = test.t;
end

C = cfg.gfx('Color_Struct');
TAG = ["X";"Y";"Z"];
fig = gobjects(2,1);

if isfield(test,'Xt')
   
   ax = gobjects(3,1);
   fig(1) = figure(...
      'Name',sprintf('%s: Preview of Cross-Trial Medians for State Data',pars.DataVariable),...
      'Color','w','Units','Normalized','Position',[1.411 0.084 0.322 0.636]);

   for ii = 1:3
      ax(ii) = subplot(3,1,ii);
      set(ax(ii),'Parent',fig(1),...
         'XColor','k','YColor','k',...
         'XGrid','on','XLim',pars.XLim,'XTick',pars.XTick,...
         'LineWidth',1.5,'FontName','Arial',...
         'NextPlot','add','Tag',TAG(ii))
      title(ax(ii),TAG(ii)); 
      ylabel(ax(ii),'median (mm)','FontName','Arial','Color','k');
      legend(ax(ii),'TextColor','black','FontName','Arial','Location','Northwest');
   end
   
   

   for ii = 1:18
      iCur = rem(ii-1,3)+1;
      plot(ax(iCur),t,test.Xt(:,ii),'Color',C.(P.Marker(ii)),...
         'LineWidth',2,'DisplayName',P.Marker(ii));
%       text(0,min(test.Xt(:,ii)),sprintf('Grand Median: %5.2f (mm)',nanmedian(test.Xt(:,ii))),...
%          'Color','k','FontName','Arial','FontWeight','bold','HorizontalAlignment','center');
   end
   xlabel(ax(3),pars.XLabel,'FontName','Arial','Color','k'); 
end


if ~isfield(test,'Xc')
   return;
end

fig(2) = figure('Name','Preview of Cross-State Medians Across Trials',...
   'Color','w','Units','Normalized','Position',[1.085 0.084 0.322 0.636]); 
nTrial = numel(unique(P.TrialID));
tAll = repmat(t,1,nTrial);
dt_jitter = min(diff(tAll))*0.01;

ax = gobjects(3,1);
for ii = 1:3
   ax(ii) = subplot(3,1,ii);
   set(ax(ii),'Parent',fig(2),...
      'XColor','k','YColor','k',...
      'XGrid','on','XLim',pars.XLim,'XTick',pars.XTick,...
      'LineWidth',1.5,'FontName','Arial',...
      'NextPlot','add','Tag',TAG(ii))
   title(ax(ii),strcat("Medoid_",TAG(ii))); 
   ylabel(ax(ii),'median (mm)','FontName','Arial','Color','k'); 
end

for ii = 1:3
   scatter(ax(ii),tAll+randn(size(tAll)).*dt_jitter,test.Xc(:,ii),...
      'MarkerFaceColor',[0.5 0.5 0.5],...
      'Marker','o',...
      'MarkerFaceAlpha',0.33,...
      'SizeData',9,...
      'MarkerEdgeAlpha',0.5);
   
   utils.addTextToAxes(ax(ii),...
      sprintf('Grand Median: %5.2f (mm)',nanmedian(test.Xc(:,ii))),...
      'southwest');
end
xlabel(ax(3),pars.XLabel,'FontName','Arial','Color','k'); 

if ~pars.AutoSave
   return;   
end

if strcmpi(pars.DataVariable,'X')
   tag = '';
else
   tag = ['-' pars.DataVariable(2:end)];
end

savefig(fig(1),fullfile('figures',[pars.CrossTrialFigureName '-' P.Properties.UserData.status tag '.fig']));
saveas(fig(1),fullfile('figures',[pars.CrossTrialFigureName '-' P.Properties.UserData.status  tag '.png']));
savefig(fig(2),fullfile('figures',[pars.CrossStateFigureName '-' P.Properties.UserData.status tag '.fig']));
saveas(fig(2),fullfile('figures',[pars.CrossStateFigureName '-' P.Properties.UserData.status tag '.png']));
delete(fig(1));
delete(fig(2));
disp('Figures saved.');

end