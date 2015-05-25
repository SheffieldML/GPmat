function [k, sk] = fileKernCompute(kern, varargin)


% FILEKERNCOMPUTE Compute the FILE kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the stored file
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG index1 : the row indices of the kernel matrix to return.
% ARG index2 : the column indices of the kernel matrix to return.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the stored file
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG index : indices of the kernel matrix to return.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : fileKernParamInit, kernCompute, kernCreate, fileKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN

sk = fileKernRead(kern, varargin{:});
k = kern.variance*sk;
