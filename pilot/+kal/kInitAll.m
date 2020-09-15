function data = kInitAll(D,k,times,N_STATE_PC)
%KINITALL Initialize multiple least-squares regression estimate and nonlinearity
%
%  Note: this is updated version of `kInit`
%
%  K = kal.kInitAll(D,k);
%  K = kal.kInitAll(D,k,times);
%
% Inputs
%  D     - Database table
%  k     - Index into row of data table
%  times - (Optional) array of times to use. If not specified, uses all.
%  N_STATE_PC - (Optional; default is 3) number of PCs to use for neural
%                 state
%
% Output
%  K     - Table for regression of each channel
%
% See also: kal, kal.showPCAreconstruction

if nargin < 4
   N_STATE_PC = 12;
end

if (nargin < 3) || isempty(times)
   times = D.Data{k}(1).times;
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

[coeffbar,Xbar,~,~,explained_bar,mubar] = pca(v);
[~,~,~,~,explained,~] = pca(Z(T,:));
coeff = coeffbar;
mu = nanmean(Z(T,:),1);
X = (Z(T,:) - mu)*coeff;
Explained_Mean = cumsum(explained_bar)./sum(explained_bar);
Explained_Trials = cumsum(explained)./sum(explained);
x = X(:,1:N_STATE_PC);
x = x-nanmean(x,1);
xbar = Xbar(:,1:N_STATE_PC);
xbar = xbar - nanmean(xbar,1);
z = Z(T,:) - nanmean(Z(T,:),1);
[A,SS_A,W] = kal.getPredictionNoise(x,xbar);
[H,SS_H,Q] = kal.getStateMeasurementNoise(z,x);
ts = repmat(tCheck,n,1);
ts = ts(T,:);
data = kal.formatData(k,0,n,T,ts,...
   coeff,mu,explained,...
   coeffbar,mubar,explained_bar,...
   N_STATE_PC,A,SS_A,W,H,SS_H,Q,Explained_Mean,Explained_Trials,xbar,x,z);


end