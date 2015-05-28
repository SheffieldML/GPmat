function [k, n2] = rbfardKernCompute(kern, x, x2)


% RBFARDKERNCOMPUTE Compute the RBFARD kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the automatic relevance determination radial basis function
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the automatic relevance determination radial basis function
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : rbfardKernParamInit, kernCompute, kernCreate, rbfardKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;
    
if nargin < 3
  n2 = dist2(x, x);
  wi2 = (.5 .* kern.inverseWidth);
  k = kern.variance*exp(-n2*wi2);
else
  x2 = x2*scales;
  n2 = dist2(x, x2);
  wi2 = (.5 .* kern.inverseWidth);
  k = kern.variance*exp(-n2*wi2);
end
