function [K_k,x_k,P,T] = getKalmanGainKinematicTimeVarying(T,varargin)
%GETKALMANGAINKINEMATICTIMEVARYING Return Kalman gain for predicting movement using spike rates, LFP, and Time-Varying Covariances
%
%  [K_k,x_k] = kal.getKalmanGainKinematicTimeVarying(T);
%
% Inputs
%  T - Main data table (init.mainData)
%  
% Output
%  K_k - Time-varying kalman gain matrix.
%  x_k - Kalman-prediction of forelimb state variables
%
% See also: kal, kal.getKalmanGainKinematic, init, init.mainData

pars = struct;
pars.CovSamplesQ = -3:3;
pars.CovSamplesW = -3:3;
pars.NIterations = 10;
pars.LFPProbe = "P1";
pars.LFPChannel = ["001","024","025"];
pars.OutputVariableName = 'Xpred_TimeVarying';
pars.OutputErrorName = 'Xerror_TimeVarying';
pars.P = [];
pars.ProcessNoiseGain = 1;
pars.ProcessInterferenceGain = 0.025;
pars.MeasurementNoiseGain = 1;
pars.ConditioningNoiseVariance = 1e-1;
fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

t = T.Properties.UserData.t;
nTrials = size(T,1);
nSamples = T.Properties.UserData.nSamples;
Tmask = cell(1,2);
[Tmask{:}] = utils.getPredictionMask(t',t',nTrials);
ich = ismember(T.Properties.UserData.Channels.Probe,pars.LFPProbe) & ...
   ismember(T.Properties.UserData.Channels.Channel,pars.LFPChannel);


Z_kin_avg = nanmean(cat(3,T.X{:}),3);
Z_all = vertcat(T.X{:}) - repmat(Z_kin_avg,nTrials,1);
Zmu = nanmean(Z_all,1);
Zsd = nanstd(Z_all,[],1);
Z_all = (Z_all - Zmu)./Zsd;


X_all = vertcat(T.Z_spike{:});
X_lfp = vertcat(T.Z_lfp{:});
X_lfp = X_lfp(:,ich);
X_xc = vertcat(T.Z_xc{:});
X_xc = X_xc(:,ich);
X_all = [X_all, X_lfp, X_xc];

X_tmp_l = cat(3,T.Z_lfp{:});
X_tmp_xc = cat(3,T.Z_xc{:});
X_xt = horzcat(cat(3,T.Z_spike{:}), X_tmp_l(:,ich,:), X_tmp_xc(:,ich,:));
Xmu_xt = nanmean(X_xt,3);
% X_all = [X_all - repmat(Xmu_xt,nTrials,1), Z_all];
X_all = X_all - repmat(Xmu_xt,nTrials,1);


% iVec = 1:nSamples;
% Xtrial_mu = cell(nTrials,1);
% for ii = 1:nTrials
%    curVec = iVec + (ii-1)*nSamples;
%    Xtrial_mu{ii} = nanmean(X_all(curVec,:),1);
%    X_all(curVec,:) = X_all(curVec,:) - Xtrial_mu{ii};
% end

Xmu = nanmean(X_all,1);
Xsd = nanstd(X_all,[],1);
X_all = (X_all - Xmu)./Xsd;

H = (Z_all(Tmask{2},:)')/(X_all(Tmask{2},:)');
A = (X_all(Tmask{2},:)')/(X_all(Tmask{1},:)');

fig = figure(...
   'Name','Generator and Relation Matrices',...
   'Color','w',...
   'Units','Normalized',...
   'Position',[0.2 0.2 0.4 0.4]);
subplot(2,1,1);
imagesc(A); 
title('A','FontName','Arial','Color','k');
colorbar;
subplot(2,1,2);
imagesc(H); 
title('H','FontName','Arial','Color','k');
colorbar;

% nLFP = sum(ich);
% nS = size(T.Z_spike{1},2);
x_k = cell(size(T,1),1);
K_k = cell(size(T,1),1);

Q = cov(Z_all(Tmask{1},:)) * pars.MeasurementNoiseGain;
W = cov(X_all(Tmask{1},:)) * pars.ProcessNoiseGain;
% Conditioning noise
% Qc = ones(size(Q)).*pars.ConditioningNoiseVariance;
% Wc = ones(size(W)).*pars.ConditioningNoiseVariance;
z_k = cell(size(T,1),1);
h = waitbar(0,'Estimating Kalman gains...');
for iTrial = 1:nTrials 
   tOffset = nSamples*(iTrial-1);
   f = (2:nSamples)+tOffset;
   b = (1:(nSamples-1))+tOffset;
   A_k = (X_all(f,:)-nanmean(X_all(f,:),1))'/(X_all(b,:)-nanmean(X_all(b,:),1))';
   H_k = (Z_all(f,:)-Z_kin_avg(2:nSamples,:)-nanmean(Z_all(f,:),1))'/(X_all(f,:)-nanmean(X_all(f,:),1))';
   Z = ((T.X{iTrial} - Z_kin_avg - Zmu)./Zsd)';
   Ko = zeros(size(H,2),size(H,1));
   P_k = eye(size(W));
%    xo = ([((T.Z_spike{iTrial}(1,:)-Xmu_xt(1,1:nS))-Xmu(1:nS))./Xsd(1:nS), ...
%           (T.Z_lfp{iTrial}(1,ich)-Xmu_xt(1,(nS+1):(nS+nLFP))-Xmu(1,(nS+1):(nS+nLFP)))./Xsd(1,(nS+1):(nS+nLFP)), ...
%           (T.Z_xc{iTrial}(1,ich)-Xmu_xt(1,(nS+nLFP+1):(nS+nLFP*2))-Xmu(1,(nS+nLFP+1):(nS+nLFP*2)))./Xsd(1,(nS+nLFP+1):(nS+nLFP*2)), ...
%           (T.X{iTrial}(1,:)-Xmu(1,(nS+nLFP*2+1):end))./Xsd(1,(nS+nLFP*2+1):end)])';
%    xo = ([((T.Z_spike{iTrial}(1,:)-Xmu_xt(1,1:nS))-Xmu(1:nS))./Xsd(1:nS), ...
%           (T.Z_lfp{iTrial}(1,ich)-Xmu_xt(1,(nS+1):(nS+nLFP))-Xmu(1,(nS+1):(nS+nLFP)))./Xsd(1,(nS+1):(nS+nLFP)), ...
%           (T.Z_xc{iTrial}(1,ich)-Xmu_xt(1,(nS+nLFP+1):(nS+nLFP*2))-Xmu(1,(nS+nLFP+1):(nS+nLFP*2)))./Xsd(1,(nS+nLFP+1):(nS+nLFP*2))])';
   xo = zeros(size(X_all,2),1);
   for iIter = 1:pars.NIterations
      K_k{iTrial} = [];
      x_k{iTrial} = [];
      Wc = zeros(size(W,1),size(W,2),size(Z,2));
      for ii = 1:size(Z,2)
%          qVec = unique(min(max(ii+pars.CovSamplesQ,1),size(Z,2)))+tOffset;
%          wVec = unique(min(max(ii+pars.CovSamplesW,1),size(Z,2)))+tOffset;
%          Q_k = cov(Z_all(qVec,:)-nanmean(Z_all(qVec,:)))*pars.MeasurementNoiseGain;
%          W_k = cov(X_all(wVec,:)-nanmean(X_all(wVec,:)))*pars.ProcessNoiseGain;
         Wc(:,:,ii) = randn(size(W)).*pars.ProcessInterferenceGain;
         W_k = W + nanmean(Wc,3);
         [Ko,xo,P_k] = kal.computeThirdOrderKf(Z(:,ii),xo,A,H,Ko,P_k,Q,W_k);
         K_k{iTrial} = cat(3,K_k{iTrial},Ko);   
         x_k{iTrial} = [x_k{iTrial}, xo];
      end
   end
  
   z_k{iTrial} = H*x_k{iTrial}.*Zsd.' + Zmu.' + Z_kin_avg.';
   x_k{iTrial} = (x_k{iTrial}.*Xsd.') + Xmu.' + Xmu_xt.';
   waitbar(iTrial/nTrials);
end
delete(h);
delete(fig);

P = pars.P;
if isempty(pars.P)
   return;
end

% z_k = cellfun(@(C)H*C,x_k,'UniformOutput',false);
P.(pars.OutputVariableName) = cell2mat(z_k);
P.(pars.OutputErrorName) = cell2mat(cellfun(@(X,Y)(X'-Y),T.X,z_k,'UniformOutput',false));
T.TimeVaryingState = x_k;
P.Properties.UserData.(pars.OutputErrorName) = T.Properties.UserData.t;
end