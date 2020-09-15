function K = kInit(D,k,times,N_STATE_PC)
%KINIT Initialize multiple least-squares regression estimate and nonlinearity
%
%  K = kal.kInit(D,k);
%  K = kal.kInit(D,k,times);
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
   N_STATE_PC = 3;
end

if nargin < 3
   times = D.Data{k}(1).times;
end

% Initialize times vectors
tCheck = D.Data{k}(1).times;
n = numel(D.Data{k});

% Get time masking vectors
T = cell(1,4);
[T{:}] = analyze.dynamics.getPredictionMask(tCheck,times,n);
ts = repmat(tCheck,n,1).*1e-3;

% Get cross-trial mean (v)
Xt = cat(3,D.Data{k}.A);
v = repmat(nanmean(cat(3,Xt),3),n,1); % Cross-trial mean
X = vertcat(D.Data{k}.A);
v_l = cell(size(T));
X_l = cell(size(T));
for iT = 1:numel(T)
   v_l{iT} = v(T{iT},:);
   X_l{iT} = X(T{iT},:) - v_l{iT};
end

% Get principal components data using cross-trial averages (remember, those
% are subtracted out prior to doing the other regression, as shown above)
Ubar = cell(1,3);
U = cell(1,3);
Explained = cell(1,3);
for ii = 1:3
   [Ubar{ii},Sbar] = svd(v(T{ii},:),'econ');
   [U{ii},~] = svd(X(T{ii},:),'econ');
   Explained{ii} = cumsum(diag(Sbar)./sum(diag(Sbar))).*100;
end


[A,SS_A,W] = kal.getPCpredictor(...
   U{1}(:,1:N_STATE_PC),U{2}(:,1:N_STATE_PC),U{3}(:,1:N_STATE_PC),...
   Ubar{1}(:,1:N_STATE_PC),Ubar{2}(:,1:N_STATE_PC),Ubar{3}(:,1:N_STATE_PC),...
   ts(T{3}));



K = [];
for ch = 1:size(X,2)
   [B,SS_B] = kal.getChannelRegression(...
      X_l{1}(:,ch),X_l{2}(:,ch),X_l{3}(:,ch),X_l{4}(:,ch),ts(T{3}));

   K = [K; ...
        table(ch,{B},... % Channel and prediction matrix
            SS_B,SS_B.Total.Rsquared,... % Fit info: H
            'VariableNames',...
            {'Channel','B','SS_B','R2_B'})]; %#ok<AGROW>
end
[B,SS_B] = kal.getChannelRegression(X_l{1},X_l{2},X_l{3},X_l{4},ts(T{3}));
[H,SS_H,Q] = kal.getPCRegression(X_l{1},X_l{2},X_l{3},...
   U{1}(:,1:N_STATE_PC),U{2}(:,1:N_STATE_PC),U{3}(:,1:N_STATE_PC),ts(T{3}));

K.Properties.UserData = struct(...
   'k',k,'nTrials',n,...
   'N_STATE_PC',N_STATE_PC,...
   'A',A,'SS_A',SS_A,'W',W,...
   'R2_A',SS_A.Total.Rsquared,...
   'B',B,'SS_B',SS_B,...
   'R2_B',SS_B.Total.Rsquared,...
   'H',H,'SS_H',SS_H,'Q',Q,...
   'R2_H',SS_H.Total.Rsquared);
K.Properties.UserData.Mask = T;
K.Properties.UserData.Explained = Explained;
end