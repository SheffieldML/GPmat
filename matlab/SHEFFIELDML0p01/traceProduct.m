function t = traceProduct(A, B)

% TRACEPRODUCT Returns the trace of the product of two matrices.
%
%	Description:
%
%	T = TRACEPRODUCT(A, B) returns the trace of the product of two
%	matrices, tr(A*B).
%	 Returns:
%	  T - the trace of the product.
%	 Arguments:
%	  A - the first matrix in the product.
%	  B - the second matrix in the product.
%	
%
%	See also
%	TRACE


%	Copyright (c) 2004 Neil D. Lawrence


t = sum(sum(A.*B'));