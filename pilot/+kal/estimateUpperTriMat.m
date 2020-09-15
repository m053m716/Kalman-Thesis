function f = estimateUpperTriMat(a,Zo,Zi,dt,Ai)
%ESTIMATEUPPERTRIMAT Function to optimize using fmincon
%
%  [f,g] = estimateUpperTriMat(a,Zo,Zi,dt,Ai)
%
% Inputs
%  a     - Vector form of regression matrix coefficients
%           -> This will be optimized
%  Zo    - Original samples used for forward-prediction
%  Zi    - "Forward-lagged" samples to be predicted
%  dt    - Timestep
%  Ai    - Matrix of indices mapping elements of a into a matrix.
%
% Output
%  f     - Objective function cost
%
% See also: kal, kal.kInit, kal.regressChannel,
%           kal.acoef2Amat, kal.Amat2acoef

A = kal.acoef2Amat(a,dt,Ai);
err = (Zi - A*Zo).^2;

f = sum(err(:));

end