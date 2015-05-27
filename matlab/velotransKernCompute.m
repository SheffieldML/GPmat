function K = velotransKernCompute(kern, varargin)

% VELOTRANSKERNCOMPUTE Compute the VELOTRANS kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the velocity translate
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel. The
% last column of the input matrix is time.
% ARG x2 : the input matrix associated with the columns of the
% kernel. The last column of the input matrix is time.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the velocity translate
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix. The last
% column of the input is time.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : velotransKernParamInit, kernCompute, kernCreate, velotransKernDiagCompute, translateKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2011

% KERN

  
for i = 1:length(varargin)
  t = varargin{i}(:, end);
  varargin{i}(:, end) = [];
  varargin{i} = varargin{i} - t*kern.velocity;
end
K = cmpndKernCompute(kern, varargin{:});
