function model = ivmOptimise(model, prior, display, innerIters, ...
			     outerIters, optimiseNoise);

% IVMOPTIMISE Optimise the IVM.

% IVM

if nargin < 6
  optimiseNoise = 1;
end
% Run IVM
for i = 1:outerIters
  model = ivmOptimiseIVM(model, display);
  model = ivmOptimiseKernel(model, prior, display, innerIters);
  if optimiseNoise
    model = ivmOptimiseIVM(model, display);
    model = ivmOptimiseNoise(model, prior, display, innerIters);
  end
  ivmDisplay(model);
end
