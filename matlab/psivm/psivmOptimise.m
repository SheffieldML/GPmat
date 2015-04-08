function models = psivmOptimise(models, prior, display, innerIters, ...
			     outerIters);

% PSIVMOPTIMISE Optimise the point set IVM.

% PSIVM

% Run IVM
for i = 1:outerIters
  models = psivmOptimiseIVM(models, display);
  models = psivmOptimiseKernel(models, prior, display, innerIters);
end
