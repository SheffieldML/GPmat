function params = biasKernExtractParam(kern)

% BIASKERNEXTRACTPARAM Extract parameters from bias kernel structure.

% IVM
if kern.linearBound
  params = log(exp(kern.variance)-1);
else
  params = log(kern.variance);
end