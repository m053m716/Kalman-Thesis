function [c0,tc] = computeZeroLagCorrelation(x,y,varargin)
%COMPUTEZEROLAGCORRELATION Compute zero-lag correlation between vector x and y
%
%  [c0,tc] = utils.computeZeroLagCorrelation(x,y,'Name',value,...);
%
% Inputs
%  x  - Column vector 1 (time-amplitude series)
%  y  - Column vector 2 (time-amplitude series)
%
%  Input vectors must be the same length.
%
% Output
%  c0 - Zero-lag expected product between the two series. 
%  tc - Time-centers of each window used to compute the expected
%        product between the two series.
%
%     Rxy(m) = E{x(n+m)*y(n)'} = E{x(n)*y(n-m)'}
%     -> Let m == 0, then this is the zero-lagged correlation.
%     -> c0 = Rxy(0), where m is set by pars.M, the window length. If
%     pars.fs is specified, then pars.M is overwritten by pars.T*pars.fs.
%     -> If pars.ts is specified, it overwrites pars.fs and pars.M (using
%        pars.T)
%
% Example
%  ```
%    test = struct('t',linspace(-0.5,0.5,100)',... % "Times"
%                  'X',[sin(linspace(0,2*pi,100))', ... % Data to correlate
%                       sin(linspace(pi/3,2*pi,100))']);
%    [test.c0,test.tc] = utils.computeZeroLagCorrelation(...
%       test.X,[],'ts',test.t);
%    figure('Name','Demo ComputeZeroLagCorrelation'); 
%    plot(test.t,test.X(:,1),'LineWidth',1.5,'DisplayName','x');
%    hold on; 
%    plot(test.t,test.X(:,2),'LineWidth',1.5,'DIsplayName','y');
%    plot(test.tc,test.c0,'LineWidth',2,'DisplayName','c0');
%  ```

pars = struct;
pars.fs = [];
pars.maxlag = 0;
pars.M = 20; % Default window length, if fs not given
pars.T = 0.2; % Default window length, if fs is gven (seconds)
pars.ts = [];
fn = fieldnames(pars);
for iV = 1:2:numel(varargin)
   idx = strcmpi(fn,varargin{iV});
   if sum(idx)==1
      pars.(fn{idx}) = varargin{iV+1};
   end
end

if isrow(x)
   x = x.';
end

if (nargin == 1) || isempty(y)
   y = x(:,end);
   x = x(:,1);
end

if numel(x)~=numel(y)
   error('Two vectors must be the same length');
end

if isrow(y)
   y = y.';
end

if ~isempty(pars.ts)
   if sum(diff(diff(pars.ts))) > eps
      error('Sample times must be monotonically increasing');
   end
   if numel(pars.ts)~=numel(x)
      error('Must provide same number of sample times as vector elements in `x` and `y`.');
   end
   pars.fs = 1/(nanmean(diff(pars.ts)));
end

if ~isempty(pars.fs)
   pars.M = round(pars.T * pars.fs);
end

N = numel(x) - pars.M + 1; % Number of "windows"
v = ((1:N) + (0:(pars.M-1))')';
% Create "lagged" matrices
X = x(v);
Y = y(v); 
c0 = nanmean(X.*Y,2); % Average dot-product for each value in window
if ~isempty(pars.ts)
   T = pars.ts(v);
   tc = nanmean(T,2);
else
   tc = nanmean(v,2);
   if ~isempty(pars.fs)
      tc = tc./pars.fs;
   end
end
end