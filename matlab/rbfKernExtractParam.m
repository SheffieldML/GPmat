function params = rbfKernExtractParam(kern)

% RBFKERNEXTRACTPARAM Extract parameters from rbf kernel structure.

% IVM

params = log([kern.inverseWidth kern.variance]);
