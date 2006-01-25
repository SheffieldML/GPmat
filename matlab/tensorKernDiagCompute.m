function k = tensorKernDiagCompute(kern, x)

% TENSORKERNDIAGCOMPUTE Compute diagonal of tensor kernel.

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
