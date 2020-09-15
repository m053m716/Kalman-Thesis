% +KAL  Package for all code applied to population analyses using Kalman filter approach
% MATLAB Version 9.7 (R2019b Update 5) 29-Aug-2020
%
% Functions
%  computeThirdOrderKf        - Compute third-older kalman filter matrix
%  estimateKF                 - Estimate kalman filter data
%  getCleanKinematics         - Apply pre-processing to lagged kinematics to clean errors due to DLC
%  getKalmanGainLFP           - Return Kalman gain for LFP 
%  getKalmanGainLFPX          - Return Kalman gain for Spikes, with LFP as state variable
%  getKalmanGainLFPxVar       - Return Kalman gain for Spikes, with LFP prediction error as state variable and forcing LTI by doing cross-trial mean-subtraction
%  getKalmanGainSpikes        - Return Kalman gain for Spikes only
%  getKalmanGainTimeVarying   - Return Kalman gain for predicting movement using spike rates, LFP, and Time-Varying Covariances
%  getPredictionNoise         - Returns matrix for using PC prior to predict observed PC
%  getStateMeasurementNoise   - Returns relationship between "states" and measurements
