function k = diagKernDiagCompute(kern, x)

% DIAGKERNDIAGCOMPUTE Compute diagonal of DIAG kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the diagonal noise covariance function kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : diagKernParamInit, kernDiagCompute, kernCreate, diagKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2011

% KERN

  if size(x, 2)>1
    error('Diag kernel requires 1-dimensional input.')
  end
  trans = str2func([kern.trans, 'Transform']);
  k = kern.variance*trans(x, 'atox');
end
