function params = linKernExtractParam(kern)

% LINKERNEXTRACTPARAM Extract parameters from linear kernel structure.

% IVM

params = log([kern.variance]);
