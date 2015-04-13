function model = fmvuParamInit(model)

% FMVUPARAMINIT FMVU model parameter initialisation.
% FORMAT
% DESC initialises the fast maximum variance unfolding
%  model structure with some default parameters.
% ARG model : the model structure which requires initialisation.
% RETURN model : the model structure with the default parameters placed in.
%
% SEEALSO : fmvuCreate, modelCreate, modelParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

[model.indices, model.D2] = findNeighbours(model.Y, model.k);
model.kappa = ones(model.N, model.k);
model.X = randn(model.N, model.q);
%model = spectralUpdateLaplacian(model);
%model = spectralUpdateX(model);
% Force updates.
param = fmvuExtractParam(model);
model = fmvuExpandParam(model, param);
