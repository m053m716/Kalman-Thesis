function F = estimateKF(data,ID)
%ESTIMATEKF Estimate kalman filter data
%
%  F = kal.estimateKF(data);
%  F = kal.estimateKF(data,ID);
%
% Inputs
%  data - Data struct from one of the kal.kInit__ functions
%  ID - (Optional) scalar numeric identifier. If not given, default is nan
%
% Output
%  F - Kalman filter data table for the experiment described by K
%
% See also: kal, kal.kInit

if nargin < 2
   ID = nan;
end

N_ITER = 3;

k = data.k;
n = data.nTrials; %
Z = data.Z'; % Corresponding neural rate signal
X = data.X'; % Corresponding neural "state" signal
s = data.samplesPerTrial; % # of samples

Q = data.Q; % Process noise covariance estimate
W = data.W; % Prediction noise estimate
H = data.H; % Measurement-state relation matrix
A = data.A; % Prediction matrix

% Format data into trials
dataZ = mat2cell(Z,size(Z,1),ones(1,n).*s);
dataX = mat2cell(X,size(X,1),ones(1,n).*s);
ts = data.ts;

F = [];
for iF = 1:n
   % Guesses for process covariance update matrix and kalman gain matrix   
   Wt = cov((dataX{iF}(:,2:end)-A*dataX{iF}(:,(1:(end-1))))');
%    Qt = cov((dataZ{iF}-H*dataX{iF})');
   Qt = Q;
   Pg = W + randn(size(W)).*1e-4;
   Kg = randn(size(X,1),size(Z,1));
   for iN = 1:N_ITER
      x_e_all = [];
      z_e_all = [];
      iCur = 0;
      for ii = 1:(s - 1)
         iCur = iCur + 1;
         [K_k,x_k,P_k] = kal.computeThirdOrderKf(dataZ{iF}(:,ii),dataX{iF}(:,ii),...
            A,H,Kg,Pg,Qt,Wt);
         x_e = dataX{iF}(:,ii+1)-x_k;
         z_e = dataZ{iF}(:,ii+1)-H*x_k;
         F = [F; ...
            table(k,iF,iN,ii,ts(iCur),{K_k},{Wt},{Qt},{x_k},{x_e},{z_e},{P_k},'VariableNames',...
               {'Block','Trial','Iteration','Sample','ts','K','Wt','Qt','xhat','xtilde','ztilde','P'})]; %#ok<AGROW>
         Kg = K_k;
         Pg = P_k;
         x_e_all = [x_e_all, x_e]; %#ok<AGROW>
         z_e_all = [z_e_all, z_e]; %#ok<AGROW>
      end
   end
end
F.Lag = repmat(data.lag,size(F,1),1);
F.ID = repmat(ID,size(F,1),1);
F = movevars(F,{'Lag','ID'},'before','ts');
F.Properties.UserData = struct('A',A,'H',H,'W',W,'Q',Q);

end