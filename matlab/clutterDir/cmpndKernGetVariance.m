function var = cmpndKernGetVariance(kern)

% CMPNDKERNGETVARIANCE Get variance of compound kernel.

var = 0;
for i = 1:length(kern.comp)
  var = var + kernGetVariance(kern.comp{i});
end
