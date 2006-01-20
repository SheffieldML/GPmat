function [Ainv, UC] = pdinv(A, UC)

% PDINV Invert a positive definite matrix.

% NDLUTIL

if nargin < 2
  UC=[];
end

% Obtain a Cholesky decomposition.
if isempty(UC)
  UC = jitChol(A);
end

invU = UC\eye(size(A, 1));
%invU = eye(size(A, 1))/UC;
Ainv = invU*invU'; 
