function K = translateKernCompute(kern, varargin)

% TRANSLATEKERNCOMPUTE Compute the TRANSLATE kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the input space translation
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the input space translation
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : translateKernParamInit, kernCompute, kernCreate,
% cmpndKernCompute, translateKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN

for i = 1:length(varargin)
  varargin{i} = varargin{i} - repmat(kern.centre, size(varargin{i}, 1), 1);
end
K = cmpndKernCompute(kern, varargin{:});
