function params = biasKernExtractParam(kern)

% BIASKERNEXTRACTPARAM Extract parameters from bias kernel structure.

% IVM

params = log(kern.variance);