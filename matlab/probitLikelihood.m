function L = probitLikelihood(X, Y, model)

% PROBITLIKELIHOOD Likelihood of data under probit noise model.

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
L = cumGaussian((Y.*mu)./(sqrt(1+varsigma)));
