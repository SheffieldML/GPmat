function params = ardKernExtractParam(kern)

% ARDKERNEXTRACTPARAM Extract parameters from ard kernel structure.

% KERN


params = [kern.inverseWidth kern.rbfVariance ...
          kern.biasVariance kern.whiteVariance ...
          kern.linearVariance kern.inputScales];
