function params = rbfKernExtractParam(kern)

% RBFKERNEXTRACTPARAM Extract parameters from rbf kernel structure.

% IVM

if kern.linearBound
  params = log(exp([kern.inverseWidth kern.variance])-1);
else
  params = log([kern.inverseWidth kern.variance]);
end