function params = ardKernExtractParam(kern)

% ARDKERNEXTRACTPARAM Extract parameters from ard kernel structure.

% IVM

if kern.linearBound
  % Optimise in the half space given by log(exp() -1)
  params = [log(exp([kern.inverseWidth kern.rbfVariance kern.biasVariance ...
                     kern.whiteVariance kern.linearVariance])-1) ...
            invSigmoid(kern.inputScales)];
else
  % Optimise in the log-space.
  params = [log([kern.inverseWidth kern.rbfVariance kern.biasVariance ...
                     kern.whiteVariance kern.linearVariance]) ...
            invSigmoid(kern.inputScales)];
end
