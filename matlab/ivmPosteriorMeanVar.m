function [mu, varsigma] = ivmPosteriorMeanVar(X, model);

% IVMPOSTERIORMEANVAR Mean and variances of the posterior at points given b X.

% IVM

D = size(model.y, 2);
numData = size(X, 1);
kX = kernCompute(X, model.kern, model.X(model.I, :))';

diagK = kernDiagCompute(X, model.kern);

varsigma = zeros(numData, D);
mu = zeros(numData, D);

if strcmp(model.noise.type, 'gaussian')
  Lk = model.Sigma.Linv*kX;
  Kinvk = model.Sigma.Linv'*Lk;
  for n = 1:numData
    varsigma(n, :) = repmat(diagK(n) - Lk(:, n)'*Lk(:, n), 1, D);
  end
end
for i = 1:D
  if ~strcmp(model.noise.type, 'gaussian')
    Lk = model.Sigma(i).Linv*kX;
    Kinvk = model.Sigma(i).Linv'*Lk;
    for n = 1:numData
      varsigma(n, i) = diagK(n) - Lk(:, n)'*Lk(:, n); 
    end
  end
  mu(:, i) = Kinvk'*model.m(model.I, i);
end