function Xout = fgplvmOptimiseSequence(model, X, Y, display, iters);

% FGPLVMOPTIMISESEQUENCE Optimise the postion of a latent sequence.
%
%	Description:
%
%	X = FGPLVMOPTIMISESEQUENCE(MODEL, X, Y, DISPLAY, ITERS) optimises
%	the location of a sequence in latent space given an initialisation
%	and an observed sequence in data space.
%	 Returns:
%	  X - the optimised sequence in the latent space.
%	 Arguments:
%	  MODEL - the model for which the point will be optimised.
%	  X - the initialisation of the sequence in the latent space.
%	  Y - the observed sequence in data space to be optimised.
%	  DISPLAY - whether or not to display the iterations of the
%	   optimisation (default: true)
%	  ITERS - maximum number of iterations for the optimisation (default
%	   2000).
%	
%
%	See also
%	FGPLVMCREATE, FGPLVMOPTIMISEPOINT, FGPLVMSEQUENCEOBJECTIVE, FGPLVMSEQUENCEGRADIENT


%	Copyright (c) 2006 Neil D. Lawrence


if nargin < 5
  iters = 2000;
  if nargin < 4
    display = true;
  end
end

options = optOptions;
if display
  options(1) = 1;
  if length(X(:))<21
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
  Xout = optim('fgplvmSequenceObjectiveGradient', X(:)', options, ...
               model, Y);
else
  % NETLAB style optimization.
  Xout = optim('fgplvmSequenceObjective', X(:)',  options, ...
                 'fgplvmSequenceGradient', model, Y);
end

Xout = reshape(Xout, size(X, 1), size(X, 2));

