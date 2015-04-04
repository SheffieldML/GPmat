function [mu, varsigma] = gpPosteriorMeanVar(model, X);

% GPPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.

% GP

D = size(model.m, 2);
numData = size(X, 1);

varsigma = zeros(numData, D);
mu = zeros(numData, D);

% Compute kernel for new point.
kX = kernCompute(model.kern, X, model.X)';

% COmpute diagonal of kernel for new point.
diagK = kernDiagCompute(model.kern, X);

Lk = model.Sigma.Linv*kX;
Kinvk = model.Sigma.Linv'*Lk;
for n = 1:numData
  varsigma(n, :) = repmat(diagK(n) - Lk(:, n)'*Lk(:, n), 1, D);
  if varsigma(n, 1) < 0
    warning('Varsigma less than zero');
  end
end
mu = Kinvk'*model.m;
