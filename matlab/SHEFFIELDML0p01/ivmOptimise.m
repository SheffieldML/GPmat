function model = ivmOptimise(model, options);

% IVMOPTIMISE Optimise the IVM.
%
%	Description:
%
%	IVMOPTIMISE(MODEL, OPTIONS) optimises an IVM by iterating between
%	selecting points and optimising kernel and noise parameters.
%	 Arguments:
%	  MODEL - the model to be optimised.
%	  OPTIONS - options structure as returned by ivmOptions.
%	
%
%	See also
%	IVMOPTIMISEIVM, IVMOPTIMISEKERNEL, IVMOPTIMISENOISE


%	Copyright (c) 2005, 2006 Neil D. Lawrence


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
