function model = ivmOptimiseNoise(model, display, iters);

% IVMOPTIMISENOISE Optimise the noise parameters.
% FORMAT
% DESC optimises the noise parameters in the IVM model.
% ARG model : the model for which the noise parameters are to be
% optimised. 
% ARG display : how much to display during optimisation (defaults
% to level 1).
% ARG iters : how many iterations of optimisation to use (defaults
% to 500).
%
% SEEALSO : optimiseParams, defaultOptions, ivmNegLogLikelihood, ivmNegGradientNoise
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

model = optimiseParams('noise', 'scg', 'ivmNegLogLikelihood', ...
                       'ivmNegGradientNoise', options, model);
