function params = linardKernExtractParam(kern)

% LINARDKERNEXTRACTPARAM Extract parameters from linear ARD kernel structure.

% IVM

params = [kern.variance kern.inputScales];
