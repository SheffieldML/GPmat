function [k, n2] = rbfKernCompute(kern, x, x2)

% RBFKERNCOMPUTE Compute the kernel given the parameters and X.
% FORMAT
% DESC computes a kernel matrix for the kernel type 'rbf' given an
% DESC input data matrix.
% ARG kern : kernel structure to be computed.
% ARG X : input data matrix (rows are data points) to the kernel computation.
% RETURN K : computed elements of the kernel structure.
%
% FORMAT 
% DESC computes a kernel matrix for the kernel type 'rbf' given an
% DESC two input data matrix, one for the rows and one for the columns.
% ARG kern : kernel structure to be computed.
% ARG X : first input matrix to the kernel computation (forms the rows of the kernel).
% ARG X2 : second input matrix to the kernel computation (forms the columns of the kernel).
% RETURN K : computed elements of the kernel structure.
%
% FORMAT
% DESC additionally returns a matrix of inter point distances for
% DESC the matrix X..
% ARG kern : kernel structure to be computed.
% ARG X : input data matrix (rows are data points) to the kernel computation.
% RETURN K : computed elements of the kernel structure.
% RETURN DIST : squared distances between the input data points.
%
% FORMAT 
% DESC additionally returns a matrix of inter point distances for
% DESC between the matrices X and X2.
% ARG kern : kernel structure to be computed.
% ARG X : first input matrix to the kernel computation (forms the rows of the kernel).
% ARG X2 : second input matrix to the kernel computation (forms the columns of the kernel).
% RETURN K : computed elements of the kernel structure.
% RETURN DIST : squared distances between the input data points.
%
% SEEALSO : kernCompute, rbfKernCreate, rbfKernDiagCompute, dist2

% KERN

if nargin < 3
  n2 = dist2(x, x);
  wi2 = (.5 .* kern.inverseWidth);
  k = kern.variance*exp(-n2*wi2);
else
  n2 = dist2(x, x2);
  wi2 = (.5 .* kern.inverseWidth);
  k = kern.variance*exp(-n2*wi2);
end
