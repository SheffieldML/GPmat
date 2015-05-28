function model = gpReversibleDynamicsExpandParam(model, params)

% GPREVERSIBLEDYNAMICSEXPANDPARAM Place the parameters vector into the model for GP dynamics.
%
%	Description:
%	model = gpReversibleDynamicsExpandParam(model, params)
%

model = gpDynamicsExpandParam(model, params);
