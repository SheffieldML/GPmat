function model = ivmOptimiseKernel(model, display, iters);

% IVMOPTIMISEKERNEL Optimise the kernel parameters.
% FORMAT
% DESC optimises the kernel parameters of the IVM model.
% ARG model : the model for which the kernel parameters are to be
% optimised. 
% ARG display : how much to display during optimisation (defaults
% to level 1).
% ARG iters : how many iterations of optimisation to use (defaults
% to 500).
%
% SEEALSO : optimiseParams, defaultOptions, ivmKernelObjective, ivmKernelGradient
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

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
  
