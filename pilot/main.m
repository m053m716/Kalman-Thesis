%MAIN Outline for repository workflow
%
% NOTE: this is moved into "pilot" because it was used initially to get the
%        Kalman formulation working. Now should use data that is all
%        relevant to the same experiment (so, LFP and Spike data as well as
%        the kinematic movement that was missing previously, which can all
%        be tied together -- this dataset is missing the kinematic
%        component which is pretty important)

clearvars -except D snips

if exist('D','var')==0
   D = getfield(load(fullfile('D:\MATLAB\Data\RC','Multi_jPCA_Table_Long-Timescale.mat'),'D'),'D');
end

if exist('snips','var')==0
   % Does not correspond to animals in dataset D: for simulation purposes
   % only! 
   % Note that pre and post is offset by 30-ms (prior) to spike data
   snips = cpl.getTrialWaveforms('E_PRE',-0.500,'E_POST',0.500);
end

% Slow oscillations are on P1-024 and P1-025 around the behavior
lfp_data = snips.SG(strcmpi(snips.Probe,'P1') & ...
   (strcmpi(snips.Channel,'024') | strcmpi(snips.Channel,'025')));
lfp_t = snips.Properties.UserData.t;

%% Organize data for kalman filter extraction
k = 116;          % Block index
ts = -300:5:300;  % Time vector (ms)
N_PC = 3;         % Number of PCs to use for neural state
LFP_LAG = -0.100:0.005:0.050;
N_ITER = 20;

% Take random subset of trial LFP data and align to these spikes using same
% alignment to behavior.
F = [];
nTotalIterations = numel(LFP_LAG)*N_ITER;
curPct = 0;
iCurrent = 0;
fprintf(1,'Synthesizing data...000%%\n');
for iLag = 1:numel(LFP_LAG)
   for iIter = 1:N_ITER
      iCurrent = iCurrent + 1;
      data = kal.kInitLFP(D,k,ts,N_PC,lfp_data,lfp_t,LFP_LAG(iLag));
      F = [F; kal.estimateKF(data,iIter)]; %#ok<AGROW>
      % F.Iteration refers to iteratively guessing the best per-trial Kalman
      % gain starting guess
      thisPct = round(iCurrent/nTotalIterations*100);
      if (thisPct-5) > curPct
         curPct = thisPct;
         fprintf(1,'\b\b\b\b\b%03d%%\n',curPct);
      end
   end
end

[G,SynthData] = findgroups(F(:,{'ID','Lag'}));
SynthData.RMS_x = splitapply(@(x)nanmean(vertcat(x{:}).^2),F.xtilde,G);
SynthData.RMS_z = splitapply(@(x)nanmean(vertcat(x{:}).^2),F.ztilde,G);

[G,TID] = findgroups(SynthData(:,'Lag'));
TID.RMS_x_mu = splitapply(@(x)nanmean(x),SynthData.RMS_x,G);
TID.RMS_x_sem = splitapply(@(x)nanstd(x)./sqrt(numel(x)),SynthData.RMS_x,G);
TID.RMS_z_mu = splitapply(@(x)nanmean(x),SynthData.RMS_z,G);
TID.RMS_z_sem = splitapply(@(x)nanstd(x)./sqrt(numel(x)),SynthData.RMS_z,G);


fig = figure('Name','Lagged LFP Error Characteristics',...
   'Color','w','Units','Normalized','Position',[0.3 0.3 0.5 0.5]);

ax = subplot(2,1,1);
set(ax,'XColor','k','YColor','k','LineWidth',1.5,...
   'FontName','Arial','NextPlot','add','Parent',fig,...
   'ColorOrder',[1 0.2 0.2; 0.6 0 0]);
bar(ax,TID.Lag,TID.RMS_x_mu);
title(ax,'\bfRMS(\itx_{tilde}\rm\bf)','FontName','Arial','Color','k');
errorbar(ax,TID.Lag,TID.RMS_x_mu,TID.RMS_x_sem,...
   'LineStyle','none','LineWidth',1.5,'Color','k');

ax = subplot(2,1,2);
set(ax,'XColor','k','YColor','k','LineWidth',1.5,...
   'FontName','Arial','NextPlot','add','Parent',fig,...
   'ColorOrder',[1 0.2 0.2; 0.6 0 0]);
bar(ax,TID.Lag,TID.RMS_z_mu);
title(ax,'\bfRMS(\itz_{tilde}\rm\bf)','FontName','Arial','Color','k');
errorbar(ax,TID.Lag,TID.RMS_z_mu,TID.RMS_z_sem,...
   'LineStyle','none','LineWidth',1.5,'Color','k');

%%

% Plot Neural State errors
fig = figure('Name','Distribution of Errors: Kalman Model Predictions',...
   'Color','w','Units','Normalized','Position',[0.2 0.2 0.6 0.6]);
nRow = floor(sqrt(N_PC+2));
nCol = ceil((N_PC+2)/nRow);
T_INTEREST = -25:5:25;
CH = ["P1-024 (LFP-synth)", "P1-025 (LFP-synth)"];
for ii = 1:(N_PC+2)
   ax = subplot(nRow,nCol,ii);
   set(ax,'XColor','k','YColor','k','LineWidth',1.5,...
      'NextPlot','add','FontName','Arial');
   ksdensity(ax,test_dist(ii,:),'Function','pdf');
   set(get(gca,'Children'),'LineWidth',2,'Color','k','Tag','T_ALL');
   ksdensity(ax,test_dist(ii,ismember(t_test,T_INTEREST)),'Function','pdf');
   set(findobj(get(gca,'Children'),'-depth',0,'Tag',''),...
      'LineWidth',1.5,'Color','b','Tag','T_INTEREST');
   if ii <= N_PC
      title(ax,sprintf('PC-%02d',ii),'FontName','Arial','Color','k');
   else
      title(ax,CH(ii-N_PC),'FontName','Arial','Color','k');
   end
   xlabel(ax,'X_{tilde}','FontName','Arial','Color','k');
   ylabel(ax,'PDF','FontName','Arial','Color','k');
end
suptitle('T: [-25 ms to +25ms]');
savefig(fig,'LFP State-Augmented Test Distributions.fig');
saveas(fig,'LFP State-Augmented Test Distributions.png');
save('LFP-State-Test.mat','test_dist','F','data','-v7.3');
delete(fig);

% Plot Neural Measurement Prediction errors
fig = figure('Name','Distribution of Measurement Prediction Errors: Kalman Model',...
   'Color','w','Units','Normalized','Position',[0.2 0.2 0.6 0.6]);
nRow = floor(sqrt(size(test_ztilde,1)));
nCol = ceil((size(test_ztilde,1))/nRow);
T_INTEREST = -25:5:25;
for ii = 1:size(test_ztilde,1)
   ax = subplot(nRow,nCol,ii);
   set(ax,'XColor','k','YColor','k','LineWidth',1.5,'NextPlot','add','FontName','Arial');
   ksdensity(ax,test_ztilde(ii,:),'Function','pdf');
   set(get(gca,'Children'),'LineWidth',2,'Color','k','Tag','T_ALL');
   ksdensity(ax,test_ztilde(ii,ismember(t_test,T_INTEREST)),'Function','pdf');
   set(findobj(get(gca,'Children'),'-depth',0,'Tag',''),'LineWidth',1.5,'Color','b','Tag','T_INTEREST');
   title(ax,sprintf('Channel-%02d',ii),'FontName','Arial','Color','k');
   xlabel(ax,'Z_{tilde}','FontName','Arial','Color','k');
   ylabel(ax,'PDF');
end
suptitle('T: [-25 ms to +25ms]');
savefig(fig,'LFP State-Augmented Measurement Prediction Error Distributions.fig');
saveas(fig,'LFP State-Augmented Measurement Prediction Error Distributions.png');
save('LFP-State-MPE_test.mat','test_ztilde','F','data','-v7.3');
delete(fig);

%% Repeat but for original data
k = 116;          % Block index
ts = -300:5:300;  % Time vector (ms)
N_PC = 3;        % Number of PCs to use for neural state
T_INTEREST = -25:5:25;
CH = ["P1-024 (LFP-synth)", "P1-025 (LFP-synth)"];

data_obs = kal.kInitAll(D,k,ts,N_PC);
F_obs = kal.estimateKF(data_obs);
obs_dist = horzcat(F_obs.xtilde{F_obs.Iteration==3});
obs_ztilde = horzcat(F_obs.ztilde{F_obs.Iteration==3});
obs_t = repmat(data_obs.ts(1:(end-1))',1,data_obs.nTrials);
fig = figure('Name','Distribution of Errors: Kalman Model Predictions',...
   'Color','w','Units','Normalized','Position',[0.2 0.2 0.6 0.6]);
nRow = floor(sqrt(N_PC));
nCol = ceil((N_PC)/nRow);
for ii = 1:(N_PC)
   ax = subplot(nRow,nCol,ii);
   set(ax,'XColor','k','YColor','k','LineWidth',1.5,'NextPlot','add','FontName','Arial');
   ksdensity(ax,obs_dist(ii,:),'Function','pdf');
   set(get(gca,'Children'),'LineWidth',2,'Color','k','Tag','T_ALL');
   ksdensity(ax,obs_dist(ii,ismember(obs_t,T_INTEREST)),'Function','pdf');
   set(findobj(get(gca,'Children'),'-depth',0,'Tag',''),'LineWidth',1.5,'Color','b','Tag','T_INTEREST');
   title(ax,sprintf('PC-%02d',ii),'FontName','Arial','Color','k');
   xlabel(ax,'X_{tilde}','FontName','Arial','Color','k');
   ylabel(ax,'PDF');
end
% savefig(fig,'Kalman Standard Prediction Error Distributions.fig');
% saveas(fig,'Kalman Standard Prediction Error Distributions.png');
% save('Standard-State-Test.mat','obs_dist','F_obs','data_obs','obs_t','-v7.3');
pause(1);
delete(fig);

fig = figure('Name','Distribution of Measurement Prediction Errors: Kalman Model',...
   'Color','w','Units','Normalized','Position',[0.2 0.2 0.6 0.6]);
nRow = floor(sqrt(size(obs_ztilde,1)));
nCol = ceil((size(obs_ztilde,1))/nRow);
T_INTEREST = -25:5:25;
for ii = 1:size(obs_ztilde,1)
   ax = subplot(nRow,nCol,ii);
   set(ax,'XColor','k','YColor','k','LineWidth',1.5,'NextPlot','add','FontName','Arial');
   ksdensity(ax,obs_ztilde(ii,:),'Function','pdf');
   set(get(gca,'Children'),'LineWidth',2,'Color','k','Tag','T_ALL');
   ksdensity(ax,obs_ztilde(ii,ismember(obs_t,T_INTEREST)),'Function','pdf');
   set(findobj(get(gca,'Children'),'-depth',0,'Tag',''),'LineWidth',1.5,'Color','b','Tag','T_INTEREST');
   title(ax,sprintf('Channel-%02d',ii),'FontName','Arial','Color','k');
   xlabel(ax,'Z_{tilde}','FontName','Arial','Color','k');
   ylabel(ax,'PDF');
end
suptitle('T: [-25 ms to +25ms]');
% savefig(fig,'Kalman Standard Measurement Prediction Error Distributions.fig');
% saveas(fig,'Kalman Standard Measurement Prediction Error Distributions.png');
% save('Standard-MPE_test.mat','obs_ztilde','F_obs','data_obs','obs_t','-v7.3');
pause(1);
delete(fig);