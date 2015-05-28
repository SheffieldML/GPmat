function H = centeringMatrix(dim)

% CENTERINGMATRIX returns the centering matrix for the given dimensionality.
% FORMAT
% DESC returns the centering matrix for the given dimensionality.
% ARG dim : the dimensionality of the centering matrix.
% RETURN H : the centering matrix for the given dimensionality.
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : eye
  
% NDLUTIL

H = -repmat(1./dim, dim, dim) + speye(dim);
