function [a,index,aeq] = Amat2acoef(A,dt,Ai)
%AMAT2ACOEF Convert from matrix to vector form for optimizer
%
%  a = kal.Amat2acoef(A,dt,Ai);
%  [a,index,aeq] = kal.Amat2acoef(A,dt,Ai);
%
% Inputs
%  A  - Matrix form of prediction matrix for channel rates
%  dt - Timestep term
%  Ai - Indexing matrix that specifies order of time-term, dt
%
% Output
%  a     - Vector form of coefficients for matrix A
%  index - Indices such that A(index) = a
%  aeq   - Equivalent indices for elements of a
%
% See also: kal, kal.estimateUpperTriMat, kal.regressChannel

n = size(A,2);
a = nan(n*(n+1)/2,1);


if nargout > 1
   index = [];
   icur = 1;
   aeq = nan(size(a));
   for ii = 1:n
      idx = find(Ai(:)==ii);
      k = numel(idx)+icur-1;
      a(icur:k) = A(idx)./(ii*dt.^ii);
      index = [index; idx]; %#ok<AGROW>
      aeq(icur:k) = ii;
      icur = k+1;
   end
else
   icur = 1;
   for ii = 1:n
      idx = find(Ai(:)==ii);
      k = numel(idx)+icur-1;
      a(icur:k) = A(idx)./(ii*dt.^ii);
      icur = k+1;
   end
end

end