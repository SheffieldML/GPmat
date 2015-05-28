function [k, sk] = whitehKernCompute(kern, x, x2)

% WHITEHKERNCOMPUTE Compute the WHITEH kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the whiteh noise
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the inpute matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the whiteh noise
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : whitehKernParamInit, kernCompute, kernCreate, whitehKernDiagCompute
%
% COPYRIGHT : Mauricio A. Alvarez, Neil D. Lawrence, 2009

% KERN

if nargin < 3
  % /~ MAURICIO : This is intended for the school Data    
  % ~/
  sk = sparseDiag(1./x(:,end));  
  k = kern.variance*sk;
else
  k = spalloc(size(x, 1), size(x2, 1), 0);
end
