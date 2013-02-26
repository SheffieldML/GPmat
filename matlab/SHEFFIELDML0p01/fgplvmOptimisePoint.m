function x = fgplvmOptimisePoint(model, x, y, display, iters);

% FGPLVMOPTIMISEPOINT Optimise the postion of a latent point.
%
%	Description:
%
%	X = FGPLVMOPTIMISEPOINT(MODEL, X, Y, DISPLAY, ITERS) optimises the
%	location of a single point in latent space given an initialisation
%	and an observed data point. Note that it ignores any dynamics
%	present in the model.
%	 Returns:
%	  X - the optimised location in the latent space.
%	 Arguments:
%	  MODEL - the model for which the point will be optimised.
%	  X - the initialisation of the point in the latent space.
%	  Y - the observed data point for which the latent point is to be
%	   optimised.
%	  DISPLAY - whether or not to display the iterations of the
%	   optimisation (default: true)
%	  ITERS - maximum number of iterations for the optimisation (default
%	   2000).
%	
%
%	See also
%	FGPLVMCREATE, FGPLVMOPTIMISESEQUENCE, FGPLVMPOINTOBJECTIVE, FGPLVMPOINTGRADIENT


%	Copyright (c) 2005, 2006 Neil D. Lawrence


if nargin < 5
  iters = 2000;
  if nargin < 4
    display = true;
  end
end

options = optOptions;
if display
  options(1) = 1;
  options(9) = 1;
end
options(14) = iters;


if isfield(model, 'optimiser')
  optim = str2func(model.optimiser);
else
  optim = str2func('scg');
end


if strcmp(func2str(optim), 'optimiMinimize')
  % Carl Rasmussen's minimize function 
  x = optim('fgplvmPointObjectiveGradient', x, options, model, y);
else
  % NETLAB style optimization.
  x = optim('fgplvmPointObjective', x,  options, ...
            'fgplvmPointGradient', model, y);
end


