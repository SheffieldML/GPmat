function g = gpGradient(params, model)

% GPGRADIENT Gradient wrapper for a GP model.
%
%	Description:
%
%	G = GPGRADIENT(PARAMS, MODEL) wraps the log likelihood gradient
%	function to return the gradient of the negative of the log
%	likelihood. This can then be used in, for example, NETLAB,
%	minimisation tools.
%	 Returns:
%	  G - the returned gradient of the negative log likelihood for the
%	   given parameters.
%	 Arguments:
%	  PARAMS - the parameters of the model.
%	  MODEL - the model for which gradients will be computed.
%	
%
%	See also
%	SCG, CONJGRAD, GPCREATE, GPGRADIENT, GPLOGLIKEGRADIENT, GPOPTIMISE


%	Copyright (c) 2005, 2006 Neil D. Lawrence


model = gpExpandParam(model, params);
g = - gpLogLikeGradients(model);
