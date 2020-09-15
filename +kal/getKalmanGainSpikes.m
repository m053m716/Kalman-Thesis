function [K_k,x_k,P] = getKalmanGainSpikes(T,varargin)
%GETKALMANGAINSPIKES Return Kalman gain for Spikes only
%
%  [K_k,x_k] = kal.getKalmanGainSpikes(T);
%
% Inputs
%  T - Main data table (init.mainData)
%  
% Output
%  K_k - Time-varying kalman gain matrix.
%  x_k - Kalman-prediction of forelimb state variables

pars = struct;
pars.NIterations = 10;
pars.P = [];
pars.OutputVariableName = 'Xpred_spikes';
pars.OutputErrorName = 'Xerror_spikes';
fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

t = T.Properties.UserData.t;
% nSamples = numel(t);
nTrials = size(T,1);
[T0,T1,T2,T3] = utils.getPredictionMask(t',t',nTrials);


Z_all = vertcat(T.Z_spike{:});
Zmu = nanmean(Z_all,1);
Zsd = nanstd(Z_all,[],1);
Z_all = (Z_all - Zmu)./Zsd;

X_all = vertcat(T.X{:});
Xmu = nanmean(X_all,1);
Xsd = nanstd(X_all,[],1);
X_all = (X_all - Xmu)./Xsd;
X_all_f = [X_all(T1,:), X_all(T2,:), X_all(T3,:)];
X_all_b = [X_all(T0,:), X_all(T1,:), X_all(T2,:)]; 

H = (Z_all(T1,:)')/(X_all_b');
Q = cov(Z_all);

A = (X_all_f')/(X_all_b');
W = cov(X_all_b);

fig = figure;
subplot(2,1,1);
imagesc(A); title('A');
subplot(2,1,2);
imagesc(H); title('H');

x_k = cell(size(T,1),1);
K_k = cell(size(T,1),1);
h = waitbar(0,'Estimating Kalman gains...');
for iTrial = 1:nTrials   
   Z = ((T.Z_spike{iTrial}-Zmu)./Zsd)';
   Ko = zeros(size(H,2),size(H,1));
   P_k = eye(size(W));
     
   xo = ([(T.X{iTrial}(1,:)-Xmu)./Xsd, (T.X{iTrial}(2,:)-Xmu)./Xsd, (T.X{iTrial}(3,:)-Xmu)./Xsd])';
   
   for iIter = 1:pars.NIterations
      K_k{iTrial} = [];
      x_k{iTrial} = [];
      for ii = 1:size(Z,2)
         [Ko,xo,P_k] = kal.computeThirdOrderKf(Z(:,ii),xo,A,H,Ko,P_k,Q,W);
         K_k{iTrial} = cat(3,K_k{iTrial},Ko);   
         x_k{iTrial} = [x_k{iTrial}, xo];
      end
   end
   
   x_k{iTrial} = x_k{iTrial}.*(repmat(Xsd,1,3)') + repmat(Xmu,1,3)';
   waitbar(iTrial/nTrials);
end
delete(h);
delete(fig);


P = pars.P;
if isempty(P)
   return;
end

tmp = cellfun(@(C)C(1:size(X_all,2),:)',x_k,'UniformOutput',false);
P.(pars.OutputVariableName) = cell2mat(cellfun(@(C)C',tmp,'UniformOutput',false));
P.(pars.OutputErrorName) = cell2mat(cellfun(@(X,Y)(X-Y)',T.X,tmp,'UniformOutput',false));

end