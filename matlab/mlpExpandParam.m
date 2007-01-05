function model = mlpExpandParam(model, params)

% MLPEXPANDPARAM Update mlp model with new vector of parameters.
% FORMAT
% DESC takes a vector of MLP weights and places them in their
% respective positions in the MLP model. The function is a wrapper
% for the mlpunpak command.
% ARG model : the model in which the weights are to be placed.
% ARG params : a vector of the weights to be placed in the model.
% RETURN model : the model with the weights distributed in the
% correct places.
%
% SEEALSO : mlpunpak, mlpCreate, mlpExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

model = mlpunpak(model, params);