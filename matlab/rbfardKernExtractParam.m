function params = rbfardKernExtractParam(kern)

% RBFARDKERNEXTRACTPARAM Extract parameters from radial basis function ARD kernel structure.

% KERN

% KERN


params = [kern.inverseWidth kern.variance kern.inputScales];
