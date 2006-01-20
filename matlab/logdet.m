function [ld, UC] = logdet(A, UC)

% LOGDET The log of the determinant when argument is positive definite.

% NDLUTIL

if nargin < 2
  UC=[];
end

% Obtain a Cholesky decomposition.
if isempty(UC)
  UC = jitChol(A);
end

ld = 2*sum(log(diag(UC)));
