function model = ivmOptimiseKernel(model, display, iters);

% IVMOPTIMISEKERNEL Optimise the kernel parameters.

% IVM

if nargin < 3
  iters = 500;
  if nargin < 2
    display = 1;
  end
end
options = defaultOptions;
if display
  options(1) = 1;
end
options(14) = iters;


model = optimiseParams('kern', 'scg', 'ivmKernelObjective', ...
                       'ivmKernelGradient', options, model);
  