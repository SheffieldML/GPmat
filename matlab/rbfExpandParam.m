function model = rbfExpandParam(model, params)

% RBFEXPANDPARAM Update rbf model with new vector of parameters.

% MLTOOLS

model = rbfunpak(model, params);