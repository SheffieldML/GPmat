function params = ardKernExtractParam(kern)

% ARDKERNEXTRACTPARAM Extract parameters from ard kernel structure.

% IVM

params = log([kern.inverseWidth kern.rbfVariance kern.biasVariance ...
              kern.whiteVariance kern.linearVariance kern.inputScales]);
