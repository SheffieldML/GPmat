function model = gpReversibleDynamicsExpandParam(model, params)

% GPREVERSIBLEDYNAMICSEXPANDPARAM Place the parameters vector into the model for GP dynamics.

% FGPLVM

model = gpDynamicsExpandParam(model, params);
