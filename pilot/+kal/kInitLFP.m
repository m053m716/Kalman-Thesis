function data = kInitLFP(D,k,times,N_STATE_PC,lfp_data,lfp_t,lfp_lag)
%KINITLFP Initialize multiple least-squares regression for LFP predictions
%
%  Note: this is updated version of `kInit`
%
%  K = kal.kInitLFP(D,k);
%  K = kal.kInitLFP(D,k,times,N_STATE_PC,lfp_data,lfp_t,lfp_lag);
%
% Inputs
%  D     - Database table
%  k     - Index into row of data table
%  times - (Optional) array of times to use. If not specified, uses all.
%  N_STATE_PC - (Optional; default is 3) number of PCs to use for neural
%                 state
%  lfp_data - LFP data, with each channel as a different cell in an array
%     and each row within cell arrays as a different trials with columns
%     representing new samples that occur at "relative" time denoted by lfp_t
%  lfp_t - Times for LFP observations
%  lfp_lag - Offset to "lag" LFP relative to rest of behavior  (SHOULD BE
%              IN SECONDS); negative lags cause the observed LFP values to
%              become shifted "earlier" relative to the occurrence of the
%              other data samples.
%
% Output
%  K     - Table for regression of each channel
%
% See also: kal, kal.showPCAreconstruction

if nargin < 7
   lfp_lag = 0;
end

if (nargin < 4) || isempty(N_STATE_PC)
   N_STATE_PC = 3;
end

if (nargin < 3) || isempty(times)
   times = D.Data{k}(1).times;
end

if (nargin < 2) || isempty(k)
   k = 116;
end

% Initialize times vectors
tCheck = D.Data{k}(1).times;
n = numel(D.Data{k});
T = repmat(ismember(tCheck,times),n,1);

% Get cross-trial mean (v)
Zt = cat(3,D.Data{k}.A);
v = repmat(nanmean(cat(3,Zt),3),n,1); % Cross-trial mean
Z = vertcat(D.Data{k}.A);

v = v(T,:);

% Get LFP at correct sample rate
ts = repmat(tCheck,n,1);
ts = ts(T,:);
ts_short = ts(1:(sum(T)/n - 1));

fs = 1/nanmean(diff(ts_short)*1e-3);
[lfp_data,lfp_t] = cpl.doDecimation(lfp_data,lfp_t+lfp_lag,round(fs));
T_LFP = round(lfp_t*1e3);
iKeep = (T_LFP >= min(ts)) & (T_LFP <= max(ts));
lfp_data = cellfun(@(C)C(:,iKeep),lfp_data,'UniformOutput',false);

nTrialSynth = size(lfp_data{1},1);
iTrial = randsample(nTrialSynth,n); % Get subset to use synthetically
X_lfp = [];
Xbar_lfp = [];
for ii = 1:numel(lfp_data)
   tmp = lfp_data{ii}(iTrial,:)';
   tmp = tmp./max(abs(tmp));
   Xbar_lfp = [Xbar_lfp, repmat(nanmean(tmp,2),n,1)]; %#ok<AGROW>
   X_lfp = [X_lfp, tmp(:)]; %#ok<AGROW>
end
[coeffbar,Xbar,~,~,explained_bar,mubar] = pca(v);
[~,~,~,~,explained,~] = pca(Z(T,:));
coeff = coeffbar;
mu = nanmean(Z(T,:),1);
X = (Z(T,:) - mu)*coeff;

Explained_Mean = cumsum(explained_bar)./sum(explained_bar);
Explained_Trials = cumsum(explained)./sum(explained);
x = [X(:,1:N_STATE_PC),X_lfp.*(1-Explained_Trials(N_STATE_PC))];
x = x - nanmean(x,1);
xbar = [Xbar(:,1:N_STATE_PC),Xbar_lfp.*(1-Explained_Mean(N_STATE_PC))];
xbar = xbar - nanmean(xbar,1);
z = Z(T,:) - nanmean(Z(T,:),1);
[A,SS_A,W] = kal.getPredictionNoise(x,xbar);
[H,SS_H,Q] = kal.getStateMeasurementNoise(z,x);
data = kal.formatData(k,lfp_lag*1e3,n,T,ts,...
   coeff,mu,explained,...
   coeffbar,mubar,explained_bar,...
   N_STATE_PC,A,SS_A,W,H,SS_H,Q,...
   Explained_Mean,Explained_Trials,xbar,x,z);

end