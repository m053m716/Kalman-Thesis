function [H,SS,Q] = getStateMeasurementNoise(Z,X)
%GETSTATEMEASUREMENTNOISE Returns relationship between "states" and measurements
%
%  [H,SS,Q] = kal.getStateMeasurementNoise(Z,X);
%
% Inputs
%  Z - Array of measured rates 
%  X - Array of estimated left singular values (neural states) for Z
%
% Output
%  H   - (Fixed / LTI) regression: "Prediction" matrix for each channel's
%         activity onto its own future values
%  SS  - Sum of squares data struct
%  Q   - Process noise covariance matrix estimate
%
% See also: kal, kal.showPCAreconstruction, kal.kInit, stat.getSS, 
%           kal.regressChannel

H = (Z')/(X');
Zhat = (H * X')';

SS = stat.getSS(Z,Zhat);
SS.x = X;
SS.npc = size(X,2);
SS.H = H;
SS.err = Z - Zhat;
Q = cov(SS.err);

end