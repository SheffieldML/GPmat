function [Ainv, UC, jitter] = pdinv(A, UC)

% PDINV Invert a positive definite matrix.
%
%	Description:
%
%	[AINV, U] = PDINV(A) inverts a positive definite matrix. If the
%	matrix isn't quite positive definite the function adds 'jitter' to
%	make it positive definite and gives out a warning message (this is
%	done through JITCHOL).
%	 Returns:
%	  AINV - the inverse of A computed using Cholesky decomposition.
%	  U - the Cholesky decomposition of A.
%	 Arguments:
%	  A - the input positive definite matrix to be inverted.
%
%	[AINV, U] = PDINV(A, U) inverts a positive definite matrix given the
%	Cholesky decomposition of A.
%	 Returns:
%	  AINV - the inverse of A computed using Cholesky decomposition.
%	  U - the Cholesky decomposition of A.
%	 Arguments:
%	  A - the input positive definite matrix to be inverted.
%	  U - the Cholesky decomposition of A.
%
%	[AINV, U, JITTER] = PDINV(A, U) inverts a positive definite matrix
%	given the Cholesky decomposition of A. If jitter is used then the
%	amount of jitter used is returned.
%	 Returns:
%	  AINV - the inverse of A computed using Cholesky decomposition.
%	  U - the Cholesky decomposition of A.
%	  JITTER - the amount of jitter added.
%	 Arguments:
%	  A - the input positive definite matrix to be inverted.
%	  U - the Cholesky decomposition of A.
%	
%
%	See also
%	JITCHOL, LOGDET, CHOL


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

invU = UC\eye(size(A, 1));
%invU = eye(size(A, 1))/UC;
Ainv = invU*invU'; 
