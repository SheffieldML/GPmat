function y = gpOut(model, x);

% GPOUT Evaluate the output of an Gaussian process model.
%
%	Description:
%
%	Y = GPOUT(MODEL, X) evaluates the output of a given Gaussian process
%	model.
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
%	GPCREATE, GPPOSTERIORMEANVAR


%	Copyright (c) 2006 Neil D. Lawrence and Carl Ek


if nargin < 2
  % This implies evaluate for the training data.
  mu = model.mu;
  varsigma = model.varSigma;
else
  if isfield(model, 'noise')
    [mu, varsigma] = gpPosteriorMeanVar(model, x);
    y = noiseOut(model.noise, mu, varsigma);
  else
    y = gpPosteriorMeanVar(model, x);
  end
end
