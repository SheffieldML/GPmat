function model = rbfExpandParam(model, params)

% RBFEXPANDPARAMS Update rbf model with new vector of parameters.

model = rbfunpak(model, params);