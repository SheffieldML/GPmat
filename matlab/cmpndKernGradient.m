function g = cmpndKernGradient(kern, x, covGrad)

% CMPNDKERNGRADIENT Gradient of compound kernel's parameters.

% IVM

g = zeros(1, kern.nParams);
startVal = 1;
endVal = 0;
for i = 1:length(kern.comp)
  endVal = endVal + kern.comp{i}.nParams;
  g(1, startVal:endVal)  = feval([kern.comp{i}.type 'KernGradient'], ...
                                      kern.comp{i}, x, covGrad);
  startVal = endVal + 1;
end
g = g*kern.paramGroups;