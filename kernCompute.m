function k = kernCompute(kern, x, x2)

% KERNCOMPUTE Compute the kernel given the parameters and X.
% FORMAT
% DESC computes a kernel matrix for the given kernel type given an
% input data matrix.
% ARG kern : kernel structure to be computed.
% ARG X : input data matrix (rows are data points) to the kernel computation.
% RETURN K : computed elements of the kernel structure.
%
% FORMAT 
% DESC computes a kernel matrix for the given kernel type given 
% two input data matrices, one for the rows and one for the columns.
% ARG kern : kernel structure to be computed.
% ARG X : first input matrix to the kernel computation (forms the rows of the kernel).
% ARG X2 : second input matrix to the kernel computation (forms the columns of the kernel).
% RETURN K : computed elements of the kernel structure.
%
% SEEALSO : kernCreate, kernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN

fhandle = str2func([kern.type 'KernCompute']);
if nargin < 3
  k = fhandle(kern, x);
else
  k = fhandle(kern, x, x2);
end
