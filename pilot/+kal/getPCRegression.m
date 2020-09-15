function [H,SS,Q] = getPCRegression(Z0,Z1,Z2,X0,X1,X2,ts)
%GETPCREGRESSION Returns per-channel regression for lagged samples onto average PC
%
%  [H,SS,Q] = kal.getPCRegression(Z0i,Z1i,Z2i,X0,X1,X2,ts);
%
%  Note: this function is redundant and is just for the sake of explicitly
%        enumerating how steps are processed. You could just pass one
%        vector with all the time-samples and lag them back and forth as
%        desired.
%
% Inputs
%  Z0i - Values at timestep (sample) 'k+0' on channel 'i', where k is an
%           index to a discrete sample, and i is an index into a spatial
%           array of microwire channels. (column vec)
%  Z1i - Values at timestep 'k+1' on channel 'i' (column vec)
%  Z2i - Values at timestep 'k+2' on channel 'i' (column vec)
%  {X0,X1,X2}   - State estimate at each timestep according to each Zi
%  ts  - Timestep at each value of X0i
%
% Output
%  H   - (Fixed / LTI) regression: "Prediction" matrix for each channel's
%         activity onto its own future values
%  SS  - Sum of squares data struct
%  Q   - Process noise covariance matrix estimate
%
% See also: kal, kal.showPCAreconstruction, kal.kInit, stat.getSS, 
%           kal.regressChannel

dt = nanmean(diff(ts));

% "Most-recent" samples at the top; lagged samples at end
Z = [Z2, (Z2-Z1)./dt, ((Z2-Z1)-(Z1-Z0))./(dt^2)];
X = [X2, (X2-X1)./dt, ((X2-X1)-(X1-X0))./(dt^2)];

H = (Z')/(X');
Z_hat = (H * X')';

SS = stat.getSS(Z,Z_hat);
SS.x = X;
SS.ts = ts;
SS.npc = size(X2,2);
SS.H = H;
SS.err = ((SS.y')-SS.H*(SS.x'))';
Q = cov(SS.err);

end