function params = cmpndKernExtractParam(kern)

% CMPNDKERNEXTRACTPARAM Extract parameters from compound kernel structure.

% IVM

params = zeros(1, kern.nParams);
startVal = 1;
endVal = 0;
for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  params(1, startVal:endVal)  = kernExtractParam(kern.comp{i});
  startVal = endVal + 1;
end
