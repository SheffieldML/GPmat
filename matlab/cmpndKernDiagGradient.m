function g = cmpndKernDiagGradient(kern, x, covDiag)

% CMPNDKERNDIAGGRADIENT Gradient of compound kernel's parameters.

% KERN

g = zeros(1, kern.nParams);
startVal = 1;
endVal = 0;

for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  if ~isempty(kern.comp{i}.index)
    % only part of the data is involved in the kernel.
    g(1, startVal:endVal) = kernDiagGradient(kern.comp{i}, ...
                                             x(:, kern.comp{i}.index), ...
                                             covDiag);
  else
    % all the data is involved with the kernel.
    g(1, startVal:endVal)  = kernDiagGradient(kern.comp{i}, x, covDiag);
  end
  startVal = endVal + 1;
end
g = g*kern.paramGroups;