function fig = singleTrialError(L,P,trialID,varargin)
%SINGLETRIALERROR Preview error from a single trial with LFP on subplot
%
%  fig = preview.singleTrialError(L,P,trialID,'Name',value,...)
%
% Inputs
%  L - LFP table as returned by L = init.lfpData();
%  P - Position table as returned by P = init.kinData(), but after it has
%     had "error" variables added to it via "getKalmanGain__" functions in
%     kal package.
%  trialID - Numeric index of trial to plot (optional)
%
% Output
%  fig - Figure handle
%
% See also: kal, preview, init.kinData, init.lfpData, kal.getKalmanGainLFP,
%           kal.getKalmanGainSpikes

pars = struct;
pars.AutoSave = false;
pars.C = cfg.gfx('Color_Struct');
pars.Var = 'spikes';
pars.PredVar = 'Xpred_%s';
pars.ErrorVar = 'Xerror_%s';
pars.LFPProbe = "P1";
pars.LFPChannel = ["024", "025"];
pars.MotorLFPProbe = "P1";
pars.MotorLFPChannel = "001";
pars.SingleTrialFigureName = 'Trial-%02d Error and LFP for %s';
pars.TROI = [-0.24 0.24];
pars.XLim = [-0.24 0.24];
pars.XLabel = 'Time (sec)';
pars.XTick = [-0.25 -0.05 0 0.05 0.25];
pars.YAxisErr = [0 35];
pars.YLimKinematic = [-20 120];
fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

if nargin < 3
   trialID = P.TrialID(randi(size(P.TrialID,1),1));
end

predVar = sprintf(pars.PredVar,pars.Var);
errVar = sprintf(pars.ErrorVar,pars.Var);

[Gp,TIDp] = findgroups(P(:,{'Marker','Dimension','Outcome'}));
mu_xt = cell2mat(splitapply(@(X){nanmean(X,1)},P.Xpred,Gp));
mu_xt = mu_xt(TIDp.Outcome==1,:);

P = P(P.TrialID==trialID,:);
iSens = ismember(L.Probe,pars.LFPProbe) & ...
        ismember(L.Channel,pars.LFPChannel);
iMotor = ismember(L.Probe,pars.MotorLFPProbe) & ...
         ismember(L.Channel,pars.MotorLFPChannel);  

L_s = L(L.TrialID==trialID & iSens,:);
L_m = L(L.TrialID==trialID & iMotor,:);  

mrk = unique(P.Marker);
lfp_s = L_s.LFP;
lfp_m = L_m.LFP;

Xp = P.(predVar)(P.Dimension=="x",:);
Yp = P.(predVar)(P.Dimension=="y",:);
Zp = P.(predVar)(P.Dimension=="z",:);

Xe = sqrt(nansum(P.(errVar)(P.Dimension=="x",:).^2,1));
Ye = sqrt(nansum(P.(errVar)(P.Dimension=="y",:).^2,1));
Ze = sqrt(nansum(P.(errVar)(P.Dimension=="z",:).^2,1));

SST = sqrt(nansum((P.Xpred - mu_xt).^2,1));
SSE = sqrt(nansum(P.(errVar).^2,1));

% m = max(abs(SST));
% SST = (SST./m)*5; % Scale to same as LFP axes
% m = max(abs(SSE));
% SSE = (SSE./m)*5; % Scale to same as LFP axes

% m = max(abs([SST,SSE]));
% SST = (SST./m)*5; % Scale to same as LFP axes
% SSE = (SSE./m)*5; % Scale to same as LFP axes

% t = P.Properties.UserData.t;
t = linspace(pars.TROI(1),pars.TROI(2),numel(SSE));

ttxt = sprintf('Trial-%02d',trialID);
fig = figure(...
   'Name',sprintf('Preview of %s Kalman Prediction Errors',ttxt),...
   'Color','w',...
   'Units','Normalized',...
   'Position',[1.32+randn(1)*0.1 0.084 0.322 0.636]); 

% For X dimension
ax = subplot(4,1,1);
yyaxis(ax,'left');
set(ax,'Parent',fig,...
      'XColor','k','YColor','k',...
      'XGrid','on','XLim',pars.XLim,'XTick',pars.XTick,...
      'LineWidth',1.5,'FontName','Arial',...
      ... 'YLim',pars.YLimKinematic,...
      'NextPlot','add');
title(ax,'Predicted_X','FontName','Arial','Color','k'); 
ylabel(ax,'X (mm)','FontName','Arial','Color','k'); 
for iMrk = 1:numel(mrk)
   h = line(ax,t,Xp(iMrk,:),...
      'Color',pars.C.(mrk(iMrk)),...
      'LineWidth',2,...
      'DisplayName',mrk(iMrk));
   h.Annotation.LegendInformation.IconDisplayStyle = 'on';
end
yyaxis(ax,'right');
ylabel(ax,'X-Error (mm)','FontName','Arial','Color','k');
h = line(ax,t,Xe,...
   'Color','r',...
   'LineWidth',1.5,...
   'LineStyle',':',...
   'DisplayName','\surd(SSE_X)');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylim(ax,pars.YAxisErr);
legend(ax,...
   'FontName','Arial',...
   'TextColor','black',...
   'EdgeColor','none',...
   'Color','w',...
   'NumColumns',4,...
   'Location','north');
utils.addTextToAxes(ax,'Markers: \it\bfLeft\rm\it Axis','northwest');
utils.addTextToAxes(ax,'Error: \it\bfRight\rm\it Axis','northeast');

% For Y-dimension
ax = subplot(4,1,2);
yyaxis(ax,'left');
set(ax,'Parent',fig,...
      'XColor','k','YColor','k',...
      'XGrid','on','XLim',pars.XLim,'XTick',pars.XTick,...
      'LineWidth',1.5,'FontName','Arial',...
      ... 'YLim',pars.YLimKinematic,...
      'NextPlot','add');
title(ax,'Predicted_Y','FontName','Arial','Color','k'); 
ylabel(ax,'Y (mm)','FontName','Arial','Color','k'); 
for iMrk = 1:numel(mrk)
   h = line(ax,t,Yp(iMrk,:),...
      'Color',pars.C.(mrk(iMrk)),...
      'LineWidth',2,...
      'DisplayName',mrk(iMrk));
   h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
yyaxis(ax,'right');
ylabel(ax,'Y-Error (mm)','FontName','Arial','Color','k');
h = line(ax,t,Ye,...
   'Color','r',...
   'LineWidth',1.5,...
   'LineStyle',':',...
   'DisplayName','\surd(SSE_{markers})');
h.Annotation.LegendInformation.IconDisplayStyle = 'on';
ylim(ax,pars.YAxisErr);
legend(ax,...
   'FontName','Arial',...
   'TextColor','black',...
   'EdgeColor','none',...
   'Color','w',...
   'NumColumns',4,...
   'Location','north');

% For Z-dimension
ax = subplot(4,1,3);
yyaxis(ax,'left');
set(ax,'Parent',fig,...
      'XColor','k','YColor','k',...
      'XGrid','on','XLim',pars.XLim,'XTick',pars.XTick,...
      'LineWidth',1.5,'FontName','Arial',...
      ... 'YLim',pars.YLimKinematic,...
      'NextPlot','add');
title(ax,'Predicted_Z','FontName','Arial','Color','k'); 
ylabel(ax,'Z (mm)','FontName','Arial','Color','k'); 
for iMrk = 1:numel(mrk)
   h = line(ax,t,Zp(iMrk,:),...
      'Color',pars.C.(mrk(iMrk)),...
      'LineWidth',2,...
      'DisplayName',mrk(iMrk));
   h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
yyaxis(ax,'right');
ylabel(ax,'Z-Error (mm)','FontName','Arial','Color','k');
ylim(ax,pars.YAxisErr);
h = line(ax,t,Ze,...
   'Color','r',...
   'LineWidth',1.5,...
   'LineStyle',':',...
   'DisplayName','\surd(SSE_Z)');
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend(ax,...
   'FontName','Arial',...
   'TextColor','black',...
   'EdgeColor','none',...
   'Color','none',...
   'NumColumns',4,...
   'Location','west');

% For LFP
ax = subplot(4,1,4);
yyaxis(ax,'left');
set(ax,'Parent',fig,...
      'XColor','k','YColor','k',...
      'XGrid','on','XLim',pars.XLim,'XTick',pars.XTick,...
      'LineWidth',1.5,'FontName','Arial',...
      ... 'YLim',pars.YLimZ,...
      'NextPlot','add');
title(ax,'LFP','FontName','Arial','Color','k'); 
ylabel(ax,L.Properties.VariableUnits{'LFP'},'FontName','Arial','Color','k'); 
xlabel(ax,pars.XLabel,'FontName','Arial','Color','k');
t = L.Properties.UserData.t;
for ii = 1:size(lfp_s,1)
   line(ax,t,lfp_s(ii,:),...
      'Color','b',...
      'LineWidth',2.0,...
      'DisplayName',sprintf('Sensory LFP %02d',ii));
end
for ii = 1:size(lfp_m,1)
   line(ax,t,lfp_m(ii,:),...
      'Color','m',...
      'LineWidth',1.5,...
      'DisplayName',sprintf('Motor LFP %02d',ii));
end
yyaxis(ax,'right');
ylabel(ax,'Total Error (mm)','FontName','Arial','Color','k');
ylim(ax,pars.YAxisErr);
if isfield(P.Properties.UserData,errVar)
   t = P.Properties.UserData.(errVar);
else
   t = P.Properties.UserData.t;
end
line(ax,t,SSE,...
   'Color','r',...
   'LineWidth',1.5,...
   'LineStyle',':',...
   'DisplayName','\surd(SSE)');
t = P.Properties.UserData.t;
line(ax,t,SST,...
   'Color','k',...
   'LineWidth',1.5,...
   'LineStyle','--',...
   'DisplayName','\surd(SST)');
legend(ax,...
   'FontName','Arial',...
   'TextColor','black',...
   'EdgeColor','none',...
   'Color','white',...
   'Location','southwest');


if strcmpi(pars.Var,'spikeslfp')
   str = 'Spikes + LFP';
elseif strcmpi(pars.Var,'spikes')
   str = 'Spikes';
else
   str = pars.Var;
end

if P.Outcome(1)==1
   str = [str, ' Kalman Errors' newline '\color{blue}\bf(Successful)'];
else
   str = [str, ' Kalman Errors' newline '\color{red}\bf(Unsuccessful)'];
end
suptitle(str);

if ~pars.AutoSave
   return;   
end

pause(0.5);

if exist('figures/Trials','dir')==0
   mkdir('figures/Trials');
end

str = sprintf(pars.SingleTrialFigureName,trialID,pars.Var);
savefig(fig,fullfile('figures','Trials',[str '.fig']));
saveas(fig,fullfile('figures','Trials',[str '.png']));
delete(fig);
disp('Figure saved.');

end