function F = estimateKF(Z,X,A,H,Q,W,varargin)
%ESTIMATEKF Estimate kalman filter data
%
%  F = kal.estimateKF(Z,X,A,H,Q,W);
%  F = kal.estimateKF(Z,X,A,H,Q,W,'Name',value,...);
%
%  Data arguments use rows to represent time-samples and columns to
%  represent variables.
%
% Inputs
%  Z     - Measurements (single-trial) matrix
%  X     - State prior estimates (single-trial) matrix
%  A     - Prediction/generator matrix
%  H     - State-Measurement relation matrix
%  Q     - Measurement (observation) noise covariance matrix
%  W     - Process noise covariance matrix
%
%  varargin: ('Name',value) pairs
%     ID - (Optional) scalar numeric identifier. Default is nan
%
% Output
%  F - Kalman filter data table for the experiment described by K
%
% See also: kal, kal.kInit


pars = struct;
pars.ID = nan;
pars.N_ITER = 3;

fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx) == 1
      pars.(fn{idx}) = varargin{iV+1};
   end
end


F = [];
% Guesses for process covariance update matrix and kalman gain matrix   
P_k = W + randn(size(W)).*1e-4;
K_k = randn(size(X,1),size(Z,1));
x_k = X(:,1);
for iN = 1:N_ITER
   x_e_all = [];
   z_e_all = [];
   iCur = 0;
   for ii = 1:(s - 1)
      iCur = iCur + 1;
      [K_k,x_k,P_k] = kal.computeThirdOrderKf(Z(:,ii),x_k,...
         A,H,K_k,P_k,Q,W);
      x_e = X(:,ii+1)-x_k;
      z_e = Z(:,ii+1)-H*x_k;
      F = [F; ...
         table(iN,ii,{K_k},{Wt},{Qt},{x_k},{x_e},{z_e},{P_k},'VariableNames',...
            {'Iteration','Sample','K','Wt','Qt','xhat','xtilde','ztilde','P'})]; %#ok<AGROW>
      x_e_all = [x_e_all, x_e]; %#ok<AGROW>
      z_e_all = [z_e_all, z_e]; %#ok<AGROW>
   end
end
F.ID = repmat(pars.ID,size(F,1),1);
F = movevars(F,{'ID'},'before',1);
F.Properties.UserData = struct('A',A,'H',H,'W',W,'Q',Q);

end