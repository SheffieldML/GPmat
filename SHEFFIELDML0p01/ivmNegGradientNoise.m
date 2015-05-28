function g = ivmNegGradientNoise(params, model)

% IVMNEGGRADIENTNOISE Wrapper function for calling noise param gradients.
%
%	Description:
%
%	G = IVMNEGGRADIENTNOISE(PARAMS, MODEL) is a wrapper function that
%	returns the negative gradients of the log likelihood with respect to
%	the noise parameters for optimisation in the NETLAB style.
%	 Returns:
%	  G - the negative gradients of the log likelihood with respect to
%	   the noise parameters.
%	 Arguments:
%	  PARAMS - the parameters where the gradients are to be computed.
%	  MODEL - the model for which the gradients are to be computed.
%	
%	
%
%	See also
%	NOISEGRADIENTPARAM, NOISEEXPANDPARAM


%	Copyright (c) 2004, 2005 Neil D. Lawrence


model.noise = noiseExpandParam(model.noise, params);
g = - noiseGradientParam(model.noise, model.mu, model.varSigma, model.y);
