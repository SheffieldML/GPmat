function model = ivmOptimise(model, prior, display, innerIters, ...
			     outerIters);

% IVMOPTIMISE Optimise the IVM.

% IVM

% Run IVM
for i = 1:outerIters
  model = ivmOptimiseIVM(model, display);
  model = ivmOptimiseNoise(model, prior, display, 100);
  model = ivmOptimiseKernel(model, prior, display, innerIters);
end
