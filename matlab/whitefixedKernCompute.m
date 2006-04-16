function k = whitefixedKernCompute(kern, x, x2)

% WHITEFIXEDKERNCOMPUTE Compute the white fixed noise kernel given the parameters and X.
% FORMAT
% DESC computes a kernel matrix for the white fixed kernel given an
% input data matrix.
% ARG kern : kernel structure to be computed.
% ARG X : input data matrix (rows are data points) to the kernel
% computation (not used in the computation, but there for compatability).
% RETURN K : computed elements of the kernel structure.
%
% FORMAT 
% DESC computes a kernel matrix for the given kernel type given 
% two input data matrices, one for the rows and one for the columns.
% ARG kern : kernel structure to be computed.
% ARG X : first input matrix to the kernel computation (forms the
% rows of the kernel, not used in the computation, but there for
% compatability).
% ARG X2 : second input matrix to the kernel computation (forms the
% columns of the kernel, not used in the computation, but there for
% compatability).
% RETURN K : computed elements of the kernel structure.
%
% SEEALSO : whiefixedKernDiagCompute
%
% COPYRIGHT : Nathaniel J. King, 2006

% KERN

if nargin < 3
  k = whiteKernCompute(kern, x);
else
  k = whiteKernCompute(kern, x, x2);
end
