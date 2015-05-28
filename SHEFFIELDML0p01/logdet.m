function [ld, UC] = logdet(A, UC)

% LOGDET The log of the determinant when argument is positive definite.
%
%	Description:
%
%	[D, U] = LOGDET(A) returns the log determinant of a positive
%	definite matrix. If the matrix isn't quite positive definite the
%	function adds 'jitter' to make it positive definite and gives out a
%	warning message (this is done through JITCHOL).
%	 Returns:
%	  D - the log determinant of A computed using Cholesky
%	   decomposition.
%	  U - the Cholesky decomposition of A.
%	 Arguments:
%	  A - the input positive definite matrix for which the log
%	   determinant is required.
%
%	[D, U, JITTER] = LOGDET(A, U) returns the log determinant of a
%	positive definite matrix given the Cholesky decomposition of A. If
%	jitter is used then the amount of jitter used is returned.
%	 Returns:
%	  D - the log determinant of A computed using Cholesky
%	   decomposition.
%	  U - the Cholesky decomposition of A.
%	  JITTER - the amount of jitter added.
%	 Arguments:
%	  A - the input positive definite matrix for which the log
%	   determinant is required.
%	  U - the Cholesky decomposition of A.
%	
%
%	See also
%	JITCHOL, PDINV, CHOL


%	Copyright (c) 2003, 2004, 2005, 2006 Neil D. Lawrence


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
