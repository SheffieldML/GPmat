function model = pmvuParamInit(model)

% PMVUPARAMINIT PMVU model parameter initialisation.
% FORMAT
% DESC initialises the probabilistic maximum variance unfolding
%  model structure with some default parameters.
% ARG model : the model structure which requires initialisation.
% RETURN model : the model structure with the default parameters placed in.
%
% SEEALSO : pmvuCreate, modelCreate, modelParamInit
%
% COPYRIGHT : Neil D. Lawrence 2009

% MLTOOLS

[model.indices, model.D2] = findNeighbours(model.Y, model.k);
model.kappa = ones(model.N, model.k);
model.gamma = 1e-4;
model.traceY = sum(sum(model.Y.*model.Y));
params = pmvuExtractParam(model);
model = pmvuExpandParam(model, params);