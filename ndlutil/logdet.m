function [ld, UC] = logdet(A, UC)

% LOGDET The log of the determinant when argument is positive definite.
% FORMAT
% DESC returns the log determinant of a positive definite matrix. If
% the matrix isn't quite positive definite the function adds 'jitter'
% to make it positive definite and gives out a warning message (this
% is done through JITCHOL).
% ARG A : the input positive definite matrix for which the log
% determinant is required.
% RETURN d : the log determinant of A computed using Cholesky
% decomposition.
% RETURN U : the Cholesky decomposition of A.
%
% FORMAT
% DESC returns the log determinant of a positive definite matrix given
% the Cholesky decomposition of A. If jitter is used then the
% amount of jitter used is returned. 
% ARG A : the input positive definite matrix for which the log
% determinant is required.
% ARG U : the Cholesky decomposition of A.
% RETURN d : the log determinant of A computed using Cholesky
% decomposition.
% RETURN U : the Cholesky decomposition of A.
% RETURN jitter : the amount of jitter added.
%
% SEEALSO : jitChol, pdinv, chol
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

ld = 2*sum(log(diag(UC)));
