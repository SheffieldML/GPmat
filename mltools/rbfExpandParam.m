function model = rbfExpandParam(model, params)

% RBFEXPANDPARAM Update rbf model with new vector of parameters.
% FORMAT
% DESC takes a vector of RBF weights and centres and places them in their
% respective positions in the RBF model. The function is a wrapper for the
% rbfunpak command.
% ARG model : the model in which the weights are to be placed.
% ARG params : a vector of the weights to be placed in the model.
% RETURN model : the model with the weights distributed in the
% correct places.
%
% SEEALSO : rbfunpak, rbfCreate, rbfExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2007

% MLTOOLS

model = rbfunpak(model, params);
