function y = ivmOut(model, x);

% IVMOUT Evaluate the output of an IVM model.
% FORMAT
% DESC evaluates the output of a given IVM model.
% ARG model : the model for which the output is being evaluated.
% ARG x : the input position for which the output is required.
% RETURN y : the output of the GP model. The function checks if
% there is a noise model associated with the GP, if there is, it is
% used, otherwise the mean of gpPosteriorMeanVar is returned.
%
% SEEALSO : ivmCreate, ivmPosteriorMeanVar, noiseOut
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2005

% IVM 

if nargin < 2
  % This implies evaluate for the training data.
  mu = model.mu;
  varsigma = model.varSigma;
else
  [mu, varsigma] = ivmPosteriorMeanVar(model, x);
end

y = noiseOut(model.noise, mu, varsigma);
