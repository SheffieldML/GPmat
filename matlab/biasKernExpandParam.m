function kern = biasKernExpandParam(kern, params)

% BIASKERNEXPANDPARAM Create kernel structure from bias's parameters.

% IVM
if kern.linearBound
  kern.variance = linearBound(params(1));
else
  kern.variance = expBound(params(1));
end