function [K_k,x_k,P_k] = computeThirdOrderKf(z,x,A,H,K,P,Q,W)
%COMPUTETHIRDORDERKF Compute third-older kalman filter matrix
%
%  [K_k,P_k] = kal.computeThirdOrderKf(z,x,A,H,K,P,Q,W);
%
% Inputs
%  z - Measured neural rates for time-step k, k-1, k-2
%  x - Approximate neural state for time-step k, k-1, k-2
%  A - State prediction matrix
%  H - State-Measurement relation matrix
%  K - Current weights on kalman gain matrix
%  P - Prediction matrix for prediction noise update
%  W - Prediction noise
%  Q - Process noise (measurement noise; noise between measurements and
%        "true" states)
%
% Output
%  K_k - Kalman gain at timestep k + 1
%  x_k - Predicted state for timestep k + 1
%  P_k - Prediction matrix at timestep k + 1
%  
% See also: kal, kal.kInit, kal.getPCpredictor, kal.getChannelRegression

x_k = A*x + K*(z - H*A*x);
P_k_minus = A*P*A' + W;
P_k = (eye(size(K,1),size(H,2)) - K*H) * P_k_minus;
K_k = P_k_minus*H'/(H*P_k_minus*H' + Q);
end