function [k, sk] = linard2KernCompute(kern, x, x2)


% LINARD2KERNCOMPUTE Compute the LINARD2 kernel given the parameters and X.
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
% SEEALSO : linard2KernParamInit, kernCompute, kernCreate, linard2KernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
%
% COPYRIGHT : Michalis K. Titsias, 2009

% KERN


    
if nargin < 3
  sk = x*sparse(diag(kern.inputScales))*x';
else  
  sk = x*sparse(diag(kern.inputScales))*x2';
end
k = sk;
