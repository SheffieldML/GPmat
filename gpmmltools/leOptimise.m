function model = leOptimise(model, display, iters)

% LEOPTIMISE Optimise an LE model.
% FORMAT
% DESC optimises a Laplacian eigenmaps model.
% ARG model : the model to be optimised.
% RETURN model : the optimised model.
%
% SEEALSO : leCreate, modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

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
