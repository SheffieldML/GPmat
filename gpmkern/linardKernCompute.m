function [k, sk] = linardKernCompute(kern, x, x2)


% LINARDKERNCOMPUTE Compute the LINARD kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the automatic relevance determination linear
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the automatic relevance determination linear
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : linardKernParamInit, kernCompute, kernCreate, linardKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009

% KERN


scales = sparse(diag(sqrt(kern.inputScales)));
x = x*scales;
    
if nargin < 3
  sk = x*x';
else
  x2 = x2*scales;
  sk = x*x2';
end
k = sk*kern.variance;
