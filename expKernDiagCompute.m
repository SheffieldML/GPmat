function k = expKernDiagCompute(kern, x)

% EXPKERNDIAGCOMPUTE Compute diagonal of EXP kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the
% exponentiated kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : expKernParamInit, kernDiagCompute, kernCreate, expKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

if kern.isStationary
  diagVal = kernDiagCompute(kern.argument, x(1, :));
  k = repmat(kern.variance*exp(diagVal)*(exp(diagVal) - 1), size(x, 1), 1);
else
  diagVal = kernDiagCompute(kern.argument, x);
  k = kern.variance*exp(diagVals).*(exp(diagVals)-1);
end
