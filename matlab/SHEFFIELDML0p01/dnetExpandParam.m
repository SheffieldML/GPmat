function model = dnetExpandParam(model, params)

% DNETEXPANDPARAM Update dnet model with new vector of parameters.
%
%	Description:
%
%	MODEL = DNETEXPANDPARAM(MODEL, PARAMS) takes a vector of DNET
%	weights and places them in their respective positions in the DNET
%	model. For single hidden layer neural networks the function is a
%	wrapper for the dnetunpak command.
%	 Returns:
%	  MODEL - the model with the weights distributed in the correct
%	   places.
%	 Arguments:
%	  MODEL - the model in which the weights are to be placed.
%	  PARAMS - a vector of the weights to be placed in the model.
%	
%
%	See also
%	DNETUNPAK, DNETCREATE, DNETEXTRACTPARAM


%	Copyright (c) 2006, 2007 Neil D. Lawrence


model.mapping = modelExpandParam(model.mapping, params(1:end-1));

func = str2func([model.betaTransform 'Transform']);
model.beta = func(params(end), 'atox');


% Update the output weights and biases.
[model.A, model.b] = modelGetOutputWeights(model.mapping);
% Update the basis if they are stored.
if model.basisStored
  [void, model.Phi] = dnetOut(model);
end

