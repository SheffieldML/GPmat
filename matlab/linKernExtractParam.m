function params = linKernExtractParam(kern)

% LINKERNEXTRACTPARAM Extract parameters from linear kernel structure.

% IVM

if kern.linearBound
  params = log(exp([kern.variance])-1);
else
  params = log([kern.variance]);
end