function model = ivmOptimise(model, options);

% IVMOPTIMISE Optimise the IVM.
% FORMAT
% DESC optimises an IVM by iterating between selecting points and
% optimising kernel and noise parameters.
% ARG model : the model to be optimised.
% ARG options : options structure as returned by ivmOptions.
%
% SEEALSO : ivmOptimiseIvm, ivmOptimiseKernel, ivmOptimiseNoise
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
 
% IVM

for i = 1:options.extIters
  if options.kernIters
    % Update the kernel if required.
    model = ivmOptimiseIvm(model, options.display);
    model = ivmOptimiseKernel(model, options.display, options.kernIters);
  end
  if options.noiseIters
    % Update the noise model if required.
    model = ivmOptimiseIvm(model, options.display);
    model = ivmOptimiseNoise(model, options.display, options.noiseIters);
  end
  if options.display
    ivmDisplay(model);
  end
end
