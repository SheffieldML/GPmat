function t = traceProduct(A, B)

% TRACEPRODUCT Returns the trace of the product of two matrices.
% FORMAT
% DESC returns the trace of the product of two matrices, tr(A*B).
% ARG A : the first matrix in the product.
% ARG B : the second matrix in the product.
% RETURN t : the trace of the product.
%
% SEEALSO : trace
%
% COPYRIGHT : Neil D. Lawrence, 2004

% NDLUTIL

t = sum(sum(A.*B'));
