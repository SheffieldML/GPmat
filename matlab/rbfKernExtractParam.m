function params = rbfKernExtractParam(kern)

% RBFKERNEXTRACTPARAM Extract parameters from rbf kernel structure.

% IVM

params = [kern.inverseWidth kern.variance];
