function model = ivmOptimiseNoise(model, display, iters);

% IVMOPTIMISENOISE Optimise the noise parameters.
%
%	Description:
%
%	IVMOPTIMISENOISE(MODEL, DISPLAY, ITERS) optimises the noise
%	parameters in the IVM model.
%	 Arguments:
%	  MODEL - the model for which the noise parameters are to be
%	   optimised.
%	  DISPLAY - how much to display during optimisation (defaults to
%	   level 1).
%	  ITERS - how many iterations of optimisation to use (defaults to
%	   500).
%	
%
%	See also
%	OPTIMISEPARAMS, DEFAULTOPTIONS, IVMNEGLOGLIKELIHOOD, IVMNEGGRADIENTNOISE


%	Copyright (c) 2004, 2005 Neil D. Lawrence



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
