function params = rbfardKernExtractParam(kern)

% RBFARDKERNEXTRACTPARAM Extract parameters from radial basis function ARD kernel structure.

% IVM


if kern.linearBound
  % Optimise in the half space given by log(exp() -1)
  params = [log(exp([kern.inverseWidth kern.variance])-1) ...
            invSigmoid(kern.inputScales)];
else
  % Optimise in the log-space.
  params = [log([kern.inverseWidth kern.variance]) invSigmoid(kern.inputScales)];
end
