function [K_k,x_k,P] = getKalmanGainLFPxVar(T,varargin)
%GETKALMANGAINLFPXVAR Return Kalman gain for Spikes, with LFP prediction error as state variable and forcing LTI by doing cross-trial mean-subtraction
%
%  [K_k,x_k] = kal.getKalmanGainLFPxVar(T);
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
pars.LFPProbe = "P1";
pars.LFPChannel = ["001","024","025"];
pars.OutputVariableName = 'Xpred_lti';
pars.OutputErrorName = 'Xerror_lti';
fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

t = T.Properties.UserData.t;
nSamples = numel(t);
nTrials = size(T,1);

Tmask = cell(1,10);
% Tmask = cell(1,4);
[Tmask{:}] = utils.getPredictionMask(t',t',nTrials);
ich = ismember(T.Properties.UserData.Channels.Probe,pars.LFPProbe) & ...
   ismember(T.Properties.UserData.Channels.Channel,pars.LFPChannel);

Z_all = vertcat(T.Z_spike{:});
Z_lfp = vertcat(T.Z_lfp{:});
Z_lfp = Z_lfp(:,ich);
Z_all = [Z_all, Z_lfp];

Z_tmp = cat(3,T.Z_lfp{:});
Z_xt = horzcat(cat(3,T.Z_spike{:}), Z_tmp(:,ich,:));
Zmu_xt = repmat(nanmean(Z_xt,3),nTrials,1);
Z_all = Z_all - Zmu_xt;

Zmu = nanmean(Z_all,1);
Zsd = nanstd(Z_all,[],1);
Z_all = (Z_all - Zmu)./Zsd;

% Z_kin_avg = cat(3,T.X{:});
% Z_kin_avg = nanmean(Z_kin_avg,3);
% Z_kin_avg = (Z_kin_avg - nanmean(Z_kin_avg,1))./std(Z_kin_avg,[],1);
% Z_kin_avg_rep = repmat(Z_kin_avg,nTrials,1);

% Z_all = [Z_all, Z_kin_avg_rep, Z_lfp];


X_xt = cat(3,T.X{:});
Xmu_xt = nanmean(X_xt,3);
X_all = vertcat(T.X{:}) - repmat(Xmu_xt,nTrials,1);
% X_lfp  = vertcat(T.Xe_lfp{:});
% X_lfp = (X_lfp(:,ich) - nanmean(X_lfp(:,ich),1))./std(X_lfp(:,ich),[],1);
% X_lfp = X_lfp(:,ich);
Xmu = nanmean(X_all,1);
Xsd = nanstd(X_all,[],1);
X_all = (X_all - Xmu)./Xsd;

% X_all_f = [X_all(Tmask{2},:), X_all(Tmask{3},:), X_all(Tmask{4},:), X_lfp(Tmask{10},:)];
% X_all_b = [X_all(Tmask{1},:), X_all(Tmask{2},:), X_all(Tmask{3},:), X_lfp(Tmask{9},:)]; 

% X_all_f = [X_all(Tmask{2},:), X_all(Tmask{3},:), X_all(Tmask{4},:), X_lfp(Tmask{2},:)];
% X_all_b = [X_all(Tmask{1},:), X_all(Tmask{2},:), X_all(Tmask{3},:), X_lfp(Tmask{1},:)]; 

X_all_f = [X_all(Tmask{8},:), X_all(Tmask{9},:), X_all(Tmask{10},:)];
X_all_b = [X_all(Tmask{1},:), X_all(Tmask{2},:), X_all(Tmask{3},:)]; 

H = (Z_all(Tmask{1},:)')/(X_all_b');
Q = cov(Z_all);

A = (X_all_b')/(X_all_f');
W = cov(X_all_b);

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

x_k = cell(size(T,1),1);
K_k = cell(size(T,1),1);
h = waitbar(0,'Estimating Kalman gains...');
for iTrial = 1:nTrials   
   Z = ([((T.Z_spike{iTrial}-Zmu_xt(1:nSamples,1:(end-sum(ich))))-Zmu(1:(end-sum(ich))))./Zsd(1:(end-sum(ich))), ...
         ((T.Z_lfp{iTrial}(:,ich)-Zmu_xt(1:nSamples,(end-sum(ich)+1):end)) - Zmu((end-sum(ich)+1):end))./Zsd((end-sum(ich)+1):end)])';
%    Z = ([(T.Z_spike{iTrial}-Zmu)./Zsd, Z_kin_avg, T.Z_lfp{iTrial}(:,ich)])';
   Ko = zeros(size(H,2),size(H,1));
   P_k = eye(size(W));
     
%    xo = ([(T.X{iTrial}(1,:)-Xmu)./Xsd, ...
%           (T.X{iTrial}(2,:)-Xmu)./Xsd, ...
%           (T.X{iTrial}(3,:)-Xmu)./Xsd, ...
%            T.Xe_lfp{iTrial}(10,ich)])';
     xo = ([((T.X{iTrial}(8,:)-Xmu_xt(8,:))-Xmu)./Xsd, ...
            ((T.X{iTrial}(9,:)-Xmu_xt(9,:))-Xmu)./Xsd, ...
            ((T.X{iTrial}(10,:)-Xmu_xt(10,:))-Xmu)./Xsd])';
   for iIter = 1:pars.NIterations
      K_k{iTrial} = [];
      x_k{iTrial} = [];
      for ii = 1:size(Z,2)
         [Ko,xo,P_k] = kal.computeThirdOrderKf(Z(:,ii),xo,A,H,Ko,P_k,Q,W);
         K_k{iTrial} = cat(3,K_k{iTrial},Ko);   
         x_k{iTrial} = [x_k{iTrial}, xo];
      end
   end
   
%    x_k{iTrial}(1:(end-sum(ich)),:) = x_k{iTrial}(1:(end-sum(ich)),:).*(repmat(Xsd,1,3)') + repmat(Xmu,1,3)';
   x_k{iTrial} = (x_k{iTrial}.*(repmat(Xsd,1,3)') + repmat(Xmu,1,3)')+ repmat(Xmu_xt,1,3)';
   waitbar(iTrial/nTrials);
end
delete(h);
delete(fig);


P = pars.P;
if isempty(P)
   return;
end

tmp = cellfun(@(C)C(1:(size(X_all,2)),:)',x_k,'UniformOutput',false);
P.(pars.OutputVariableName) = cell2mat(cellfun(@(C)C',tmp,'UniformOutput',false));
P.(pars.OutputErrorName) = cell2mat(cellfun(@(X,Y)(X-Y)',T.X,tmp,'UniformOutput',false));

end