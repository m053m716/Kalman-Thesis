function A = acoef2Amat(a,dt,Ai)
%ACOEF2AMAT Convert from vector to matrix form for optimizer
%
%  A = kal.acoef2Amat(a,dt,Ai);
%
% Inputs
%  a  - Vector form of coefficients for matrix A
%  dt - Timestep term
%  Ai - Indexing matrix that specifies order of time-term, dt
%
% Output
%  A  - Matrix form of prediction matrix for channel rates
%
% See also: kal, kal.estimateUpperTriMat, kal.regressChannel

A = zeros(size(Ai));
icur = 1;
for ii = 1:size(Ai,2)
   idx = Ai(:)==ii;
   k = sum(idx)+icur-1;
   A(idx) = a(icur:k).*(dt^ii);
   icur = k+1;
end

end