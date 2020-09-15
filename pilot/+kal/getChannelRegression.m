function [B,SS] = getChannelRegression(Z0,Z1,Z2,Z3,ts)
%GETCHANNELREGRESSION Returns per-channel regression for lagged samples
%
%  [B,SS] = kal.getChannelRegression(Z0i,Z1i,Z2i,Z3i,ts);
%
%  Note: this function is redundant and is just for the sake of explicitly
%        enumerating how steps are processed. You could just pass one
%        vector with all the time-samples and lag them back and forth as
%        desired.
%
% Inputs
%  Z0 - Values at timestep (sample) 'k+0' on channel 'i', where k is an
%           index to a discrete sample, and i is an index into a spatial
%           array of microwire channels. (column vec)
%  Z1 - Values at timestep 'k+1' on channel 'i' (column vec)
%  Z2 - Values at timestep 'k+2' on channel 'i' (column vec)
%  Z3 - Values at timestep 'k+3' on channel 'i' (column vec)
%  dt  - Timestep at each value
%
% Output
%  B   - (Fixed / LTI) regression: "Prediction" matrix for each channel's
%         activity onto its own future values
%
% See also: kal, kal.showPCAreconstruction, kal.kInit, stat, stat.getSS, 
%           kal.regressChannel

dt = nanmean(diff(ts));

% "Most-recent" samples at the top; lagged samples at end
Zf = [Z3, (Z3-Z2)./dt, ((Z3-Z2)-(Z2-Z1))./(dt^2)]; 
Zb = [Z2, (Z2-Z1)./dt, ((Z2-Z1)-(Z1-Z0))./(dt^2)];

B = (Zf')/(Zb');
Zf_hat = (B * Zb')';

SS = stat.getSS(Zf,Zf_hat);
SS.x = Zb;
SS.ts = ts;

end