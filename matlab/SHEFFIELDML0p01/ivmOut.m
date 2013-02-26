function y = ivmOut(model, x);

% IVMOUT Evaluate the output of an IVM model.
%
%	Description:
%
%	Y = IVMOUT(MODEL, X) evaluates the output of a given IVM model.
%	 Returns:
%	  Y - the output of the GP model. The function checks if there is a
%	   noise model associated with the GP, if there is, it is used,
%	   otherwise the mean of gpPosteriorMeanVar is returned.
%	 Arguments:
%	  MODEL - the model for which the output is being evaluated.
%	  X - the input position for which the output is required.
%	
%
%	See also
%	IVMCREATE, IVMPOSTERIORMEANVAR, NOISEOUT


%	Copyright (c) 2003, 2005 Neil D. Lawrence


if nargin < 2
  % This implies evaluate for the training data.
  mu = model.mu;
  varsigma = model.varSigma;
else
  [mu, varsigma] = ivmPosteriorMeanVar(model, x);
end

y = noiseOut(model.noise, mu, varsigma);
