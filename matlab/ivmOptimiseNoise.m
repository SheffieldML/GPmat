function model = ivmOptimiseNoise(model, prior, display, iters);

% IVMOPTIMISENOISE Optimise the noise parameters.

% IVM

if nargin < 4
  iters = 500;
  if nargin < 3
    display = 1;
    if nargin < 2
      prior = 0;
    end
  end
end
options = foptions;
if display
  options(1) = 1;
end
options(14) = iters;

model = optimiseParams('noise', 'scg', 'negNoiseLogLikelihood', ...
                       'negNoiseGradientParam', options, model, prior);
