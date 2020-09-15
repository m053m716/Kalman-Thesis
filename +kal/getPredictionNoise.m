function [A,SS,W] = getPredictionNoise(X,Xbar)
%GETPREDICTIONNOISE  Return "prior" matrix (generator for state predictions) as well as state noise covariance
%
%  [A,SS,W] = kal.getPredictionNoise(X0,X1,X2,Xbar0,Xbar1,Xbar2,ts);
%
%  Note: this function is redundant and is just for the sake of explicitly
%        enumerating how steps are processed. You could just pass one
%        vector with all the time-samples and lag them back and forth as
%        desired.
%
% Inputs
%  X           - Kinematic state
%  Xbar        - Kinematic state "prior"
%
% Output
%  A   - (Fixed / LTI) regression: "Prediction" matrix for each channel's
%         activity onto its own future values
%  SS  - Sum of squares data struct
%  W   - Prediction noise covariance matrix estimate
%
% See also: kal, kal.showPCAreconstruction, kal.kInit, stat.getSS, 
%           kal.regressChannel

A = (Xbar')/(X');
Xhat = (A * X')';

SS = stat.getSS(Xbar,Xhat);
SS.x = X;
SS.npc = size(Xbar,2);
SS.A = A;
SS.err = Xbar - Xhat;
W = cov(SS.err);

end