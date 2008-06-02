function y = gpOut(model, x);

% GPOUT Evaluate the output of an Gaussian process model.
% FORMAT
% DESC evaluates the output of a given Gaussian process model.
% ARG model : the model for which the output is being evaluated.
% ARG x : the input position for which the output is required.
% RETURN y : the output of the GP model. The function checks if
% there is a noise model associated with the GP, if there is, it is
% used, otherwise the mean of gpPosteriorMeanVar is returned.
%
% SEEALSO : gpCreate, gpPosteriorMeanVar
%
% COPYRIGHT : Neil D. Lawrence and Carl Ek, 2006

% GP 

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
