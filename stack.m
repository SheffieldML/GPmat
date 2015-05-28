function x = stack(X)

% STACK Return column stacked vector of given matrix.
% FORMAT
% DESC returns a column stacked vector of a given matrix. This is
% useful if you wish to stack a column vector from a matrix
% returned by another function (i.e. you can't apply the colon
% operator directly).
% ARG X : the matrix to be stacked.
% RETURN x : stacked column vector from the matrix.
%
% SEEALSO : colon
%
% COPYRIGHT : Neil D. Lawrence, 2004


% NDLUTIL

x = X(:);
