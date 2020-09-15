function [P,Z_k,X_k,Z_e,X_e] = getCleanKinematics(P,varargin)
%GETCLEANKINEMATICS Apply pre-processing to lagged kinematics to clean errors due to DLC
%
%  P = kal.getCleanKinematics(P,varargin);
%
% Inputs
%  P        - Table of kinematic position data
%  varargin - (Optional) 'Name',value input argument pairs for parameters
%
% Output
%  P     - 'Clean' table of kinematic position data
%
% See also: kal, kal.showPCAreconstruction, kal.kInit, stat, stat.getSS, 
%           kal.regressChannel, utils.getKinematicStates

pars = struct;
pars.MakeFigures = false; 
pars.RejectionThresholdSD = 5;
pars.TROI = [-0.240 0.240]; % "Time" region-of-interest

fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx) == 1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

[Z,Xc,Xt,P] = utils.getKinematicStates(P,'DoMedianSubtraction',false);
dt = 1/P.Properties.UserData.fs;

Zmu = nanmean(Z,1);
Z = Z - Zmu;


nSamples = size(Xt,1);
nTrials = numel(unique(P.TrialID));

tCheck = (pars.TROI(1):dt:pars.TROI(2))';
tvec = P.Properties.UserData.t';
T = cell(1,4);
[T{:}] = utils.getPredictionMask(tCheck,tvec,nTrials);
Z_l = cell(size(T));
% X_l = cell(size(T));
for ii = 1:numel(T)
   Z_l{ii} = Z(T{ii},:);
%    X_l{ii} = Xc(T{ii},:);
%    X_l{ii} = Xc(T{ii},:);
end
nSamples_roi = numel(tCheck);

% Get prediction for motion
Xf = [Z_l{2}, Z_l{3}, Z_l{4}];
Xb = [Z_l{1}, Z_l{2}, Z_l{3}];
A = (Xf')/(Xb');
Xf_hat = (A * Xb')';
SS = stat.getSS(Xf,Xf_hat);
R2_A = SS.Total.Rsquared;
Xf_tilde = Xf - Xf_hat;
W = cov(Xf_tilde); % Prediction variance: noise due to factors unrelated to Newtonian motion
% W = cov(Xf);

% (Optional): Make figures
if pars.MakeFigures
   fig = figure('Name','Generator Matrix','Color','w');
   ax = subplot(2,1,1);
   set(ax,'Parent',fig,'FontName','Arial',...
      'XColor','none','YColor','none','NextPlot','add',...
      'YDir','reverse','Color','none');
   imagesc(ax,A);
   xlabel(ax,'State-Dimension_i','FontName','Arial','Color','k');
   ylabel(ax,'State-Dimension_j','FontName','Arial','Color','k');
   colorbar(ax);
   title(ax,'\bfA\rm (\itX_f = A * X_b + W\rm | a.u.)',...
      'FontName','Arial','Color','k','FontWeight','bold');
   ax = subplot(2,1,2);
   set(ax,'Parent',fig,'FontName','Arial',...
      'XColor','none','YColor','none','NextPlot','add',...
      'YDir','reverse','Color','none');
   imagesc(ax,W);
   xlabel(ax,'State-Dimension_i','FontName','Arial','Color','k');
   ylabel(ax,'State-Dimension_j','FontName','Arial','Color','k');
   colorbar(ax);
   title(ax,'\bfW\rm (Process Noise | mm^2)',...
      'FontName','Arial','Color','k','FontWeight','bold');
   savefig(fig,'figures/Matrix A - Prediction of Movement.fig');
   saveas(fig,'figures/Matrix A - Prediction of Movement.png');
   delete(fig);
end

% Get relation between "state" of hand and the individual measurements of
% the markers.
% Zb = [Z_l{1}, Z_l{2}, Z_l{3}];
Zf = Z_l{4};
H = (Zf')/(Xb');
Z_hat = (H * Xb')';
SS = stat.getSS(Zf,Z_hat);
R2_H = SS.Total.Rsquared;
Z_tilde = Zf - Z_hat;
Q = cov(Z_tilde); % Process variance: noise on our kinematic measurements
% Q = cov(Zb);

% (Optional): Make figures
if pars.MakeFigures
   fig = figure('Name','State Relation Matrix','Color','w'); 
   ax = subplot(2,1,1);
   set(ax,'Parent',fig,'FontName','Arial',...
      'XColor','none','YColor','none','NextPlot','add',...
      'YDir','reverse','Color','none');
   imagesc(ax,H);
   xlabel(ax,'Marker-Variable_i','FontName','Arial','Color','k');
   ylabel(ax,'State-Dimension_j','FontName','Arial','Color','k');
   colorbar(ax);
   title(ax,'\bfH\rm (\itZ_{hat} = H * X_c + Q\rm | a.u.)',...
      'FontName','Arial','Color','k','FontWeight','bold');
   ax = subplot(2,1,2);
   set(ax,'Parent',fig,'FontName','Arial',...
      'XColor','none','YColor','none','NextPlot','add',...
      'YDir','reverse','Color','none');
   imagesc(ax,Q);
   xlabel(ax,'Marker-Variable_i','FontName','Arial','Color','k');
   ylabel(ax,'Marker-Variable_j','FontName','Arial','Color','k');
   colorbar(ax);
   title(ax,'\bfQ\rm (Observation Noise | mm^2)',...
      'FontName','Arial','Color','k','FontWeight','bold');
   savefig(fig,'figures/Matrix H - Hand Medoid State to Marker Relations.fig');
   saveas(fig,'figures/Matrix H - Hand Medoid State to Marker Relations.png');
   delete(fig);
end

T = cell(1,3);
[T{:}] = utils.getPredictionMask(tvec,tvec,nTrials);
% Z_l = cell(1,3);
X_l = cell(1,3);
for ii = 1:3
%    Z_l{ii} = Z(T{ii},:);
   X_l{ii} = Xc(T{ii},:); 
end

% Z_obs = [Z_l{1},Z_l{2},Z_l{3}];
Z_obs = Z(T{1},:);
% X_obs = [X_l{1},X_l{2},X_l{3}];
X_obs = [Z(T{1},:),Z(T{2},:),Z(T{3},:)];

% Final step is to actually clean the data using the kalman formulation.
% Do this for each trial, individually: each time, intialize xhat as the
% first value in our reduced-roi xhat model expectation.
Z_tr = cellfun(@(C)C',mat2cell(Z_obs,ones(1,nTrials).*(nSamples-2),size(Zf,2)),'UniformOutput',false);
X_tr = cellfun(@(C)C',mat2cell(X_obs,ones(1,nTrials).*(nSamples-2),size(Xb,2)),'UniformOutput',false);
h = waitbar(0,'Applying Kalman filter to clean kinematic data (preprocessing step)...');
% Kalman estimated outputs:
X_k = [];
Z_k = []; 
X_e = [];
Z_e = [];

for iTrial = 1:nTrials
   x_k = Xb(1,:)';
   K_k = zeros(size(H,2),size(H,1));
   P_k = eye(size(W));
   X_k = [X_k; x_k'];
   z_k = H*x_k;
   z_e = Z_tr{iTrial}(:,1) - z_k;
   x_e = X_tr{iTrial}(:,1) - x_k;
   Z_k = [Z_k; z_k'];
   Z_e = [Z_e; z_e'];
   X_e = [X_e; x_e'];
   for iTime = 1:(nSamples-2)
      [K_k,x_k,P_k] = kal.computeThirdOrderKf(Z_tr{iTrial}(:,iTime),x_k,...
         A,H,K_k,P_k,Q,W);
      z_k = H*x_k;      
      Z_k = [Z_k; z_k']; %#ok<*AGROW>
      X_k = [X_k; x_k'];
      
      if iTime < (nSamples-2)
         z_e = Z_tr{iTrial}(:,iTime+1) - z_k;
         x_e = X_tr{iTrial}(:,iTime+1) - x_k;
      else
         Z_k = [Z_k; z_k'];
         X_k = [X_k; x_k']; % Use last prediction twice
         z_e = nan(size(z_e,1),2);
         x_e = nan(size(x_e,1),2);
      end
      Z_e = [Z_e; z_e'];
      X_e = [X_e; x_e'];
      
   end
   waitbar(iTrial/nTrials);
end
delete(h);

P.Properties.UserData.nTrials = nTrials;
P.Properties.UserData.nSamples = nSamples;
P.Properties.UserData.nSamples_roi = nSamples_roi;
P.Properties.UserData.Q = Q;
P.Properties.UserData.W = W;
P.Properties.UserData.Rsquared = struct('A',R2_A,'H',R2_H);

if nargout ~= 1
   return;
end

% Check for times when the data is a certain number of deviations from the
% predicted kalman output, and in those instances, swap the observed data
% Z_kf = Z_k(:,1:(size(Z_k,2)/3)); % Only take the "current step" states.
% Z_kf = Z_kf + repmat(X_k(:,1:3),1,size(Z_kf,2)/3);
Z_kf = Z_k(:,1:size(Z_k,2)) + Zmu; % Only take the "current step" states.
Z = Z + Zmu;
sd_z = nanstd(Z,[],1);
iBadObs = false(size(Z));
for iZ = 1:size(Z,2)
   iBadObs(:,iZ) = abs(Z(:,iZ)-Z_kf(:,iZ)) > (pars.RejectionThresholdSD*sd_z(iZ));
   Z(iBadObs(:,iZ),iZ) = Z_kf(iBadObs(:,iZ),iZ);
end


% Only apply the cleaning if the only thing to be returned is the original
% table.
Zc = mat2cell(Z,ones(1,nTrials).*nSamples,size(Z,2));
Zfiltall = mat2cell(Z_kf,ones(1,nTrials).*nSamples,size(Z_kf,2));
iBadObs = mat2cell(iBadObs,ones(1,nTrials).*nSamples,size(Z,2));

% Remember, this will ultimately become our "state" variable:
Xout = cell2mat(cellfun(@(C)C',Zc,'UniformOutput',false));
XfiltAll = cell2mat(cellfun(@(C)C',Zfiltall,'UniformOutput',false));
iBad = cell2mat(cellfun(@(C)C',iBadObs,'UniformOutput',false));

P.X = Xout;
P.Xpred = XfiltAll;
iKeep = (P.Properties.UserData.t > pars.TROI(1)) & (P.Properties.UserData.t < pars.TROI(2));
P.X = P.X(:,iKeep);
P.Xpred = P.Xpred(:,iKeep);
P.ObservationState = iBad(:,iKeep);
P.Properties.UserData.status = 'clean';
P.Properties.UserData.t = P.Properties.UserData.t(iKeep);

end