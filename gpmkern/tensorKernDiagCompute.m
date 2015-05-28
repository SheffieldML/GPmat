function k = tensorKernDiagCompute(kern, x)


% TENSORKERNDIAGCOMPUTE Compute diagonal of TENSOR kernel.
% FORMAT
% DESC computes the diagonal of the kernel matrix for the tensor product kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
% RETURN k : a vector containing the diagonal of the kernel matrix
% computed at the given points.
%
% SEEALSO : tensorKernParamInit, kernDiagCompute, kernCreate, tensorKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN


i = 1;
if ~isempty(kern.comp{i}.index)
  % only part of the data is involved with the kernel.
  k  = kernDiagCompute(kern.comp{i}, x(:, kern.comp{i}.index));
else
  % all the data is involved with the kernel.
  k  = kernDiagCompute(kern.comp{i}, x);
end
for i = 2:length(kern.comp)
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved with the kernel.
    k  = k.*kernDiagCompute(kern.comp{i}, x(:, kern.comp{i}.index));
  else
    % all the data is involved with the kernel.
    k  = k.*kernDiagCompute(kern.comp{i}, x);
  end
end
