function g = modelOutputGradX(model, X)

% MODELOUTPUTGRADX Compute derivatives with respect to model inputs of model outputs.
%
%	Description:
%
%	G = MODELOUTPUTGRADX(MODEL, X) gives the gradients of the outputs
%	from the model with respect to the inputs.
%	 Returns:
%	  G - gradients of the model output with respect to the model
%	   parameters for the given input locations.
%	 Arguments:
%	  MODEL - the model structure for which gradients are computed.
%	  X - input locations where gradients are to be computed.
%	
%
%	See also
%	MODELCREATE, MODELOUTPUTGRAD, MODELLOGLIKELIHOOD, MODELLOGLIKEGRADIENTS


%	Copyright (c) 2006 Neil D. Lawrence


fhandle = str2func([model.type 'OutputGradX']);
g = fhandle(model, X);

