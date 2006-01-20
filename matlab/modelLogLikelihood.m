function model = modelLogLikelihood(model)

% MODELLOGLIKELIHOOD Compute a model log likelihood.

% MLTOOLS

fhandle = str2func([model.type 'LogLikelihood']);
model = fhandle(model);