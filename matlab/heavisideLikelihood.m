function L = heavisideLikelihood(X, Y, model)

% HEAVISIDELIKELIHOOD Likelihood of data under heaviside noise model.

% IVM
D = size(model.y, 2);
if isempty(X)
  mu = model.mu;
  varsigma = model.varSigma;
  Y = model.y;
else
  [mu, varsigma] = ivmPosteriorMeanVar(X, model);
end
for i = 1:D
  mu(:, i) = mu(:, i) + model.noise.bias(i);
end
L = (1-2*model.noise.eta)*cumGaussian((Y.*mu)./(sqrt(varsigma)))+model.noise.eta;
