function model = leOptimise(model, display, iters)

% LEOPTIMISE Optimise an LE model.
%
%	Description:
%
%	MODEL = LEOPTIMISE(MODEL) optimises a Laplacian eigenmaps model.
%	 Returns:
%	  MODEL - the optimised model.
%	 Arguments:
%	  MODEL - the model to be optimised.
%	
%
%	See also
%	LECREATE, MODELOPTIMISE


%	Copyright (c) 2009 Neil D. Lawrence


switch model.weightType
 case 'constant'
  model.kappa = repmat(1, model.N, model.k);
 case 'rbf'
  model.kappa = exp(-dist2(model.Y, model.Y)/(2*model.weightScale*model.weightScale));
 otherwise
  error('Unknown weight type in leOptimise');
end
model = spectralUpdateLaplacian(model);
model = spectralUpdateX(model);
