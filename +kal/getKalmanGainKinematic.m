function [K_k,x_k] = getKalmanGainKinematic(T,varargin)
%GETKALMANGAINKINEMATIC Return Kalman gain for predicting movement using spike rates and LFP
%
%  [K_k,x_k] = kal.getKalmanGainKinematic(T);
%
% Inputs
%  T - Main data table (init.mainData)
%  
% Output
%  K_k - Time-varying kalman gain matrix.
%  x_k - Kalman-prediction of forelimb state variables

pars = struct;
pars.NIterations = 10;
pars.LFPProbe = "P1";
pars.LFPChannel = ["001","024","025"];
pars.OutputVariableName = 'Xpred_kinematics';
pars.OutputErrorName = 'Xerror_kinematics';
pars.ProcessNoiseGain = 1;
pars.MeasurementNoiseGain = 1;
fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

t = T.Properties.UserData.t;
nTrials = size(T,1);
Tmask = cell(1,3);
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
X_all = [X_all - repmat(Xmu_xt,nTrials,1), Z_all];

Xmu = nanmean(X_all,1);
Xsd = nanstd(X_all,[],1);
X_all = (X_all - Xmu)./Xsd;


H = (Z_all(Tmask{3},:)')/(X_all(Tmask{2},:)');
Q = cov(Z_all(Tmask{3},:)) * pars.MeasurementNoiseGain;

A = (X_all(Tmask{2},:)')/(X_all(Tmask{1},:)');
% W = cov(X_all(Tmask{2},:));
W = cov(X_all(Tmask{2},:)) * pars.ProcessNoiseGain;

fig = figure(...
   'Name','Generator and Relation Matrices',...
   'Color','w',...
   'Units','Normalized',...
   'Position',[0.2 0.2 0.4 0.4]);
subplot(2,1,1);
imagesc(A); title('A');
colorbar;
subplot(2,1,2);
imagesc(H); title('H');
colorbar;

nLFP = sum(ich);
nS = size(T.Z_spike{1},2);
x_k = cell(size(T,1),1);
% z_k = cell(size(T,1),1);
K_k = cell(size(T,1),1);
h = waitbar(0,'Estimating Kalman gains...');
for iTrial = 1:nTrials   
   Z = ((T.X{iTrial} - Zmu)./Zsd)';
   Ko = zeros(size(H,2),size(H,1));
   P_k = eye(size(W));
   xo = ([((T.Z_spike{iTrial}(1,:)-Xmu_xt(1,1:nS))-Xmu(1:nS))./Xsd(1:nS), ...
          (T.Z_lfp{iTrial}(1,ich)-Xmu_xt(1,(nS+1):(nS+nLFP))-Xmu(1,(nS+1):(nS+nLFP)))./Xsd(1,(nS+1):(nS+nLFP)), ...
          (T.Z_xc{iTrial}(1,ich)-Xmu_xt(1,(nS+nLFP+1):(nS+nLFP*2))-Xmu(1,(nS+nLFP+1):(nS+nLFP*2)))./Xsd(1,(nS+nLFP+1):(nS+nLFP*2)), ...
          (T.X{iTrial}(1,:)-Xmu(1,(nS+nLFP*2+1):end))./Xsd(1,(nS+nLFP*2+1):end)])';
   for iIter = 1:pars.NIterations
      K_k{iTrial} = [];
      x_k{iTrial} = [];
      for ii = 1:size(Z,2)
         [Ko,xo,P_k] = kal.computeThirdOrderKf(Z(:,ii),xo,A,H,Ko,P_k,Q,W);
         K_k{iTrial} = cat(3,K_k{iTrial},Ko);   
         x_k{iTrial} = [x_k{iTrial}, xo];
      end
   end
  
   x_k{iTrial} = (x_k{iTrial}.*Xsd') + Xmu' + [Xmu_xt, Z_kin_avg]';
   waitbar(iTrial/nTrials);
end
delete(h);
delete(fig);

end