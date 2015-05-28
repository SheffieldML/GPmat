function g = dnetGradient(params, model)

% DNETGRADIENT Density Network gradient wrapper.
%
%	Description:
%
%	G = DNETGRADIENT(PARAMS, MODEL) is a wrapper function for the
%	gradient of the negative log likelihood of an Density Network model
%	with respect to the latent postions and parameters.
%	 Returns:
%	  G - the gradient of the negative log likelihood with respect to
%	   the latent positions and the parameters at the given point.
%	 Arguments:
%	  PARAMS - vector of parameters and latent postions where the
%	   gradient is to be evaluated.
%	  MODEL - the model structure into which the latent positions and
%	   the parameters will be placed.
%	
%
%	See also
%	DNETLOGLIKEGRADIENTS, DNETEXPANDPARAM


%	Copyright (c) 2008 Neil D. Lawrence


model = dnetExpandParam(model, params);
g = - dnetLogLikeGradients(model);
