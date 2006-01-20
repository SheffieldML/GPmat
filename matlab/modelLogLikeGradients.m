function g = modelLogLikeGradients(model)

% MODELLOGLIKELIHOODGRADIENTS Compute a model's gradients wrt log likelihood.

% MLTOOLS

fhandle = str2func([model.type 'LogLikeGradients']);
g = fhandle(model);