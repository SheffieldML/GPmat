function model = fgplvmOptimise(model, display, iters, g_check);

% FGPLVMOPTIMISE Optimise the FGPLVM.
%
%	Description:
%
%	MODEL = FGPLVMOPTIMISE(MODEL, DISPLAY, ITERS) takes a given GP-LVM
%	model structure and optimises with respect to parameters and latent
%	positions.
%	 Returns:
%	  MODEL - the optimised model.
%	 Arguments:
%	  MODEL - the model to be optimised.
%	  DISPLAY - flag dictating whether or not to display optimisation
%	   progress (set to greater than zero) (default value 1).
%	  ITERS - number of iterations to run the optimiser for (default
%	   value 2000).
%	fgplvmLogLikeGradients, fgplvmObjective, fgplvmGradient
%	
%
%	See also
%	FGPLVMCREATE, FGPLVMLOGLIKELIHOOD, 


%	Copyright (c) 2005, 2006 Neil D. Lawrence


if(nargin<4)
  g_check = false;
  if nargin < 3
    iters = 2000;
    if nargin < 2
      display = 1;
    end
  end
end

params = fgplvmExtractParam(model);

options = optOptions;
if(display)
  options(1) = 1;
  if(length(params) <= 100||g_check)
    options(9) = 1;
  end
end
options(14) = iters;

if isfield(model, 'optimiser')
  optim = str2func(model.optimiser);
else
  optim = str2func('scg');
end

if strcmp(func2str(optim), 'optimiMinimize')
  % Carl Rasmussen's minimize function 
  params = optim('fgplvmObjectiveGradient', params, options, model);
else
  % NETLAB style optimization.
  params = optim('fgplvmObjective', params,  options, ...
                 'fgplvmGradient', model);
end

model = fgplvmExpandParam(model, params);
