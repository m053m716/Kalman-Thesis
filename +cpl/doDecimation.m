function [data_dec,t_dec] = doDecimation(lfp_data,lfp_t,fs)
%DODECIMATION Decimate data to specified sample rate
%
%  [data_dec,t_dec] = cpl.doDecimation(lfp_data,lfp_t);
%  [data_dec,t_dec] = cpl.doDecimation(lfp_data,lfp_t,fs);
%
% Inputs
%  lfp_data - Cell array of LFP trial-aligned data from desired channels.
%              One cell element per channel.
%  lfp_t    - Corresponding sample times relative to alignment for each
%              column of elements in lfp_data. Rows of each cell element in
%              lfp_data are new trials.
%  fs       - (Optional) fs desired for decimation.
%
% Output
%  data_dec - Same cell structure as lfp_data, but with decimated data
%  t_dec    - Times for sample values in data_dec
%
% See also: cpl, kal, cpl.getTrialWaveforms

if nargin < 3
   fs = 50; % Sample rate to decimate to
end

if iscell(lfp_data)
   data_dec = cell(size(lfp_data));
   for ii = 1:numel(lfp_data)
      [data_dec{ii},t_dec] = cpl.doDecimation(lfp_data{ii},lfp_t,fs);
   end
   return;
end

[x,t] = resample(lfp_data',lfp_t',fs);
data_dec = x';
t_dec = t';

end