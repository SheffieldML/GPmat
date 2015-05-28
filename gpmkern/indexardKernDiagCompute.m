function k = indexardKernDiagCompute(kern, x)

% INDEXARDKERNDIAGCOMPUTE Compute diagonal of INDEXARD kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the index ard based covariance function kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : indexKernParamInit, kernDiagCompute, kernCreate, indexKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2011

% KERN
  k = zeros(size(x, 1), 1);
  for i = 1:length(kern.indices)
    ind = find(round(x)==kern.indices(i));
    k(ind) = kern.indexScales(i);
  end
end
