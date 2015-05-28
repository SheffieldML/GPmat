function k = diagKernCompute(kern, x, x2)

% DIAGKERNCOMPUTE Compute the DIAG kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the diagonal noise covariance function
% kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 : the input matrix associated with the columns of the kernel.
% RETURN k : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the diagonal noise covariance function
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : the kernel matrix computed at the given points.
%
% SEEALSO : diagKernParamInit, kernCompute, kernCreate, diagKernDiagCompute
%
% COPYRIGHT : Neil D. Lawrence, 2011

% KERN
  if size(x, 2)>1
    error('Diag kernel requires 1-dimensional input.')
  end
  trans = str2func([kern.trans, 'Transform']);
  if nargin < 3
    sk = diag(trans(x, 'atox'));
    k = kern.variance*sk;
    k = sparse(k);
  else
    k = spalloc(size(x, 1), size(x2, 1), 0);
  end
end
