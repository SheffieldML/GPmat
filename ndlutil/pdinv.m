function [Ainv, UC, jitter] = pdinv(A, UC)

% PDINV Invert a positive definite matrix.
% FORMAT
% DESC inverts a positive definite matrix. If the matrix isn't
% quite positive definite the function adds 'jitter' to make it
% positive definite and gives out a warning message (this is done
% through JITCHOL).
% ARG A : the input positive definite matrix to be inverted.
% RETURN Ainv : the inverse of A computed using Cholesky
% decomposition.
% RETURN U : the Cholesky decomposition of A.
%
% FORMAT
% DESC inverts a positive definite matrix given the Cholesky
% decomposition of A.
% ARG A : the input positive definite matrix to be inverted.
% ARG U : the Cholesky decomposition of A.
% RETURN Ainv : the inverse of A computed using Cholesky
% decomposition.
% RETURN U : the Cholesky decomposition of A.
%
% FORMAT
% DESC inverts a positive definite matrix given the Cholesky
% decomposition of A. If jitter is used then the
% amount of jitter used is returned. 
% ARG A : the input positive definite matrix to be inverted.
% ARG U : the Cholesky decomposition of A.
% RETURN Ainv : the inverse of A computed using Cholesky
% decomposition.
% RETURN U : the Cholesky decomposition of A.
% RETURN jitter : the amount of jitter added.
%
% SEEALSO : jitChol, logdet, chol
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2005, 2006


% NDLUTIL

if nargin < 2
  UC=[];
end

% Obtain a Cholesky decomposition.
if isempty(UC)
  if nargout > 2
    [UC, jitter] = jitChol(A);
  else
    UC = jitChol(A);
  end
end

invU = UC\eye(size(A, 1));
%invU = eye(size(A, 1))/UC;
Ainv = invU*invU'; 
