function [A,SS,Z,T] = regressChannel(X,mu,s,nLag,dt)
%REGRESSCHANNEL Returns per-channel regression for lagged samples
%
%  [A,SS,Z,T] = kal.regressChannel(X,mu,s,nLag,dt);
%
%  More generalized form of kal.getChannelRegression
%
% Inputs
%  X     - Single-channel values across all time samples/recordings
%  mu    - Cross-trial mean, replicated to matched samples of X
%  s     - Number of (retained) samples per trial
%  nLag  - Number of lags to consider
%  dt    - Timestep (seconds; fixed)
%
% Output
%  A   - (Fixed / LTI) regression: "Prediction" matrix for each channel's
%         activity onto its own future values
%  SS  - Struct with sum-of-squares info about fit
%  Z   - Cell array of lagged samples used in fit, starting with Z0
%  T   - Cell array of corresponding mask vectors used on full sample
%        vector X for each matched cell element of Z.
%
% See also: kal, kal.showPCAreconstruction, kal.kInit, stat.getSS,
%           kal.getChannelRegression, analyze.dynamics.getPredictionMask

n = size(X,1)/s; % Number of trials
T = cell(1,nLag);
iTimes = true(s,1);

[T{:}] = analyze.dynamics.getPredictionMask(iTimes,n);
Z = cell(size(T));
for iT = 1:numel(T)
   Z{iT} = X(T{iT},:) - mu(T{iT},:);
end

% "Most-recent" samples at the top; lagged samples at end
Zi = horzcat(Z{2:end}); 
Zo = horzcat(Z{1:(end-1)});

% A = (Zi')/(Zo');
Ai = zeros(nLag-1);
for iLag = 1:nLag
   Ai = Ai + diag(ones(nLag-iLag,1)*iLag,iLag-1);
end
[a0,idx,~] = kal.Amat2acoef(Ai,dt,Ai);
% Equivalent:
% a0 = repelem(((1:(nLag-1))./((dt.^(1:(nLag-1))).*(1:(nLag-1))))',(nLag-1):-1:1);

[iRow,~] = ind2sub(size(Ai),idx);

bcon = ones(nLag-1,1);
Acon = zeros(nLag-1,numel(a0));
icur = 1;
for iLag = 1:(nLag-1)
   Acon(iLag,iRow==iLag) = ones(1,sum(iRow==iLag));
   bcon(iLag) = sum(a0(iRow==iLag));
   icur = icur + (nLag-iLag) + 1;
end


a = fmincon(@(a)kal.estimateUpperTriMat(a,Zo',Zi',dt,Ai),a0,Acon,bcon);
A = kal.acoef2Amat(a,dt,Ai);


Zi_hat = (A * Zo')';

SS = stat.getSS(Zi,Zi_hat);
SS.x = Zo;



end