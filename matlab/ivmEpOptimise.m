function model = ivmEpOptimise(model, prior, display, innerIters, ...
			     outerIters, optimiseNoise);

% IVMEPOPTIMISE Optimise the IVM making use of point removal.


if nargin < 6
  optimiseNoise = 1;
end
% Run IVM
for i = 1:outerIters
  model = ivmSelectPoints(model, display);
  model = ivmOptimiseKernel(model, prior, display, innerIters);
  if optimiseNoise
    model = ivmSelectPoints(model, display);
    model = ivmOptimiseNoise(model, prior, display, innerIters);
  end
  ivmDisplay(model);
end
