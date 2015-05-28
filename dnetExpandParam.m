function model = dnetExpandParam(model, params)

% DNETEXPANDPARAM Update dnet model with new vector of parameters.
% FORMAT
% DESC takes a vector of DNET weights and places them in their
% respective positions in the DNET model. For single hidden layer
% neural networks the function is a wrapper for the dnetunpak command.
% ARG model : the model in which the weights are to be placed.
% ARG params : a vector of the weights to be placed in the model.
% RETURN model : the model with the weights distributed in the
% correct places.
%
% SEEALSO : dnetunpak, dnetCreate, dnetExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2007

% MLTOOLS

model.mapping = modelExpandParam(model.mapping, params(1:end-1));

func = str2func([model.betaTransform 'Transform']);
model.beta = func(params(end), 'atox');


% Update the output weights and biases.
[model.A, model.b] = modelGetOutputWeights(model.mapping);
% Update the basis if they are stored.
if model.basisStored
  [void, model.Phi] = dnetOut(model);
end

