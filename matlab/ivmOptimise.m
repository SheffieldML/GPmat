function model = ivmOptimise(model, options);

% IVMOPTIMISE Optimise the IVM.

% IVM

% Run IVM
for i = 1:options.extIters
  if options.kernIters
    % Update the kernel if required.
    model = ivmOptimiseIVM(model, options.display);
    model = ivmOptimiseKernel(model, options.display, options.kernIters);
  end
  if options.noiseIters
    % Update the noise model if required.
    model = ivmOptimiseIVM(model, options.display);
    model = ivmOptimiseNoise(model, options.display, options.noiseIters);
  end
  ivmDisplay(model);
end
