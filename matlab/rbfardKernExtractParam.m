function params = rbfardKernExtractParam(kern)

% RBFARDKERNEXTRACTPARAM Extract parameters from radial basis function ARD kernel structure.

% KERN


params = [kern.inverseWidth kern.variance kern.inputScales];
