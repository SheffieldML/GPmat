function params = linardKernExtractParam(kern)

% LINARDKERNEXTRACTPARAM Extract parameters from linear ARD kernel structure.

% IVM

if kern.linearBound
  params = [log(exp(kern.variance)-1) invSigmoid(kern.inputScales)];
else
  params = [log(kern.variance) invSigmoid(kern.inputScales)];
end