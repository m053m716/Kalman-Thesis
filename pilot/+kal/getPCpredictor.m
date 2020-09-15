function [A,SS,W] = getPCpredictor(X0,X1,X2,Xbar0,Xbar1,Xbar2,ts)
%GETPCPREDICTOR Returns matrix for using lagged PC prior to predict observed PC
%
%  [A,SS,W] = kal.getPCpredictor(X0,X1,X2,Xbar0,Xbar1,Xbar2,ts);
%
%  Note: this function is redundant and is just for the sake of explicitly
%        enumerating how steps are processed. You could just pass one
%        vector with all the time-samples and lag them back and forth as
%        desired.
%
% Inputs
%  {X0,X1,X2}            - Neural state (as estimated by left singular values) 
%  {Xbar0,Xbar1,Xbar2}   - Neural state prior: average population left singular values (U)
%  ts                    - Timestep at each value of X0i
%
% Output
%  A   - (Fixed / LTI) regression: "Prediction" matrix for each channel's
%         activity onto its own future values
%  SS  - Sum of squares data struct
%  W   - Prediction noise covariance matrix estimate
%
% See also: kal, kal.showPCAreconstruction, kal.kInit, stat.getSS, 
%           kal.regressChannel

dt = nanmean(diff(ts));

% "Most-recent" samples at the top; lagged samples at end
Xo = [X2, (X2-X1)./dt, ((X2-X1)-(X1-X0))./(dt^2)];
Xbar = [Xbar2, (Xbar2-Xbar1)./dt, ((Xbar2-Xbar1)-(Xbar1-Xbar0))./(dt^2)];

A = (Xbar')/(Xo');
X_hat = (A * Xo')';

SS = stat.getSS(Xbar,X_hat);
SS.x = Xo;
SS.ts = ts;
SS.npc = size(Xbar2,2);
SS.A = A;
SS.err = ((SS.y')-SS.A*(SS.x'))';
SS.err_c = SS.err .* [ones(1,SS.npc), ones(1,SS.npc).*dt, ones(1,SS.npc).*dt^2];
W = cov(SS.err);

end