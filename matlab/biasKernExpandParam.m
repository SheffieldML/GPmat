function kern = biasKernExpandParam(params, kern)

% BIASKERNEXPANDPARAM Create kernel structure from bias's parameters.

% IVM

kern.variance = exp(params(1));