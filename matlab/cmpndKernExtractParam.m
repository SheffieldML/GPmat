function params = cmpndKernExtractParam(kern)

% CMPNDKERNEXTRACTPARAM Extract parameters from compound kernel structure.

% KERN

params = zeros(1, kern.nParams);
startVal = 1;
endVal = 0;
for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  params(1, startVal:endVal) = kernExtractParam(kern.comp{i});
  startVal = endVal + 1;
end
paramGroups = kern.paramGroups;
for i = 1:size(paramGroups, 2)
  ind = find(paramGroups(:, i));
  paramGroups(ind(2:end), i) = 0;
end
params = params*paramGroups;
