function kern = rbfardKernExpandParam(kern, params)

% RBFARDKERNEXPANDPARAM Create kernel structure from radial basis function ARD's parameters.

% IVM
if kern.linearBound
  params(1:2) = linearBound(params(1:2));
else
  params(1:2) = expBound(params(1:2));
end
kern.inverseWidth = params(1);
kern.variance = params(2);
kern.inputScales = sigmoidBound(params(3:end));
