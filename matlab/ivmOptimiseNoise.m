function model = ivmOptimiseNoise(model, display, iters);

% IVMOPTIMISENOISE Optimise the noise parameters.

% IVM


if nargin < 3
  iters = 500;
  if nargin < 2
    display = 1;
  end
end
options = foptions;
if display
  options(1) = 1;
end
options(14) = iters;

model = optimiseParams('noise', 'scg', 'ivmNegLogLikelihood', ...
                       'ivmNegGradientNoise', options, model);
