function [K, argK, kDiagMat] = expKernCompute(kern, varargin)

% EXPKERNCOMPUTE Compute the EXP kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the exponentiated
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the exponentiated
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : expKernParamInit, kernCompute, kernCreate, expKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN


argK = kernCompute(kern.argument, varargin{:});
if isfield(kern.argument, 'isStationary') & kern.argument.isStationary
  kDiagMat = kernDiagCompute(kern.argument, varargin{1}(1, :));
  K = kern.variance*exp(kDiagMat)*(exp(argK) - 1);
else
  if length(varargin)>1
    kii = kernDiagCompute(kern.argument, varargin{1});
    kjj = kernDiagCompute(kern.argument, varargin{2});
    kDiagMat = repmat(kii', size(kjj, 1), 1) ...
        + repmat(kjj, 1, size(kii, 1));
    kDiagMat = kDiagMat/2;
  else
    kii = kernDiagCompute(kern.argument, varargin{1});
    kDiagMat = repmat(kii', size(kii, 1), 1) ...
        + repmat(kii, 1, size(kii, 1));
    kDiagMat = kDiagMat/2;
  end
  K = kern.variance*(exp(argK)-1).*exp(kDiagMat);
end
