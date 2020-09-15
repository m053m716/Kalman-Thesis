% +KAL  Package for all code applied to population analyses using Kalman filter approach
% MATLAB Version 9.7 (R2019b Update 5) 29-Aug-2020
%
%  Table references:
%     D = utils.loadTables('multi');
%     E = analyze.dynamics.exportTable(D);
%     E2 = analyze.dynamics.exportSubTable(D);
%
%  Functions (Table D)
%     acoef2Amat                 - Convert from vector to matrix form for optimizer
%     Amat2acoef                 - Convert from matrix to vector form for optimizer
%     computeThirdOrderKf        - Compute third-older kalman filter matrix
%     estimateUpperTriMat        - Function to optimize using fmincon
%     formatData                 - Put data struct into stereotyped struct format
%     getPCpredictor             - Returns matrix for using lagged PC prior to predict observed PC
%     getPredictionNoise         - Returns matrix for using PC prior to predict observed PC
%     getStateMeasurementNoise   - Returns relationship between "states" and measurements
%     kInitAll                   - Initialize multiple least-squares regression estimate and nonlinearity
%     regressChannel             - Returns per-channel regression for lagged samples
%
%  Graphics (Table D)
%     animateKF                  - Return animation stack for a given trial to export filter weights video output
%     showPCAreconstruction      - Make bar graph of channelwise error from PCA reconstruction
