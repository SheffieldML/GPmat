function params = rbfKernExtractParam(kern)

% RBFKERNEXTRACTPARAM Extract parameters from rbf kernel structure.

% KERN

% KERN


params = [kern.inverseWidth kern.variance];
