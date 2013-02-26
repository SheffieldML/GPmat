function [f, g] = gpObjectiveGradient(params, model)

% GPOBJECTIVEGRADIENT Wrapper function for GP objective and gradient.
%
%	Description:
%
%	[F, G] = GPOBJECTIVEGRADIENT(PARAMS, MODEL) returns the negative log
%	likelihood of a Gaussian process model given the model structure and
%	a vector of parameters. This allows the use of NETLAB minimisation
%	functions to find the model parameters.
%	 Returns:
%	  F - the negative log likelihood of the GP model.
%	  G - the gradient of the negative log likelihood of the GP model
%	   with respect to the parameters.
%	 Arguments:
%	  PARAMS - the parameters of the model for which the objective will
%	   be evaluated.
%	  MODEL - the model structure for which the objective will be
%	   evaluated.
%	
%
%	See also
%	MINIMIZE, GPCREATE, GPGRADIENT, GPLOGLIKELIHOOD, GPOPTIMISE


%	Copyright (c) 2005, 2006 Neil D. Lawrence


% Check how the optimiser has given the parameters
if size(params, 1) > size(params, 2)
  % As a column vector ... transpose everything.
  transpose = true;
  model = gpExpandParam(model, params');
else
  transpose = false;
  model = gpExpandParam(model, params);
end

f = - gpLogLikelihood(model);
if nargout > 1
  g = - gpLogLikeGradients(model);
end
if transpose
  g = g';
end