function k = cmpndKernDiagCompute(kern, x)

% CMPNDKERNDIAGCOMPUTE Compute diagonal of compound kernel.
% FORMAT
% DESC computes the diagonal of a kernel matrix for a compound
% kernel type given an input data matrix.
% ARG kern : kernel structure to be computed.
% ARG X : input data matrix (rows are data points) to the kernel computation.
% RETURN K : vector containing computed diagonal elements of the
% kernel structure.
%
% SEEALSO : kernDiagCompute, cmpndKernCompute

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
    k  = k + kernDiagCompute(kern.comp{i}, x(:, kern.comp{i}.index));
  else
    % all the data is involved with the kernel.
    k  = k + kernDiagCompute(kern.comp{i}, x);
  end
end
