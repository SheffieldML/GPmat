function [params, transform] = sqexpKernExtractParam(kern)

% SQEXPKERNEXTRACTPARAM Extract parameters from squared exponential kernel structure.

% IVM

if kern.linearBound
  % Optimise in the half space given by log(exp() -1)
  params = [log(exp([kern.inverseWidth kern.rbfVariance kern.biasVariance ...
                     kern.whiteVariance])-1)];
else
  % Optimise in the log-space.
  params = log([kern.inverseWidth kern.rbfVariance ...
                kern.biasVariance kern.whiteVariance]);
end

