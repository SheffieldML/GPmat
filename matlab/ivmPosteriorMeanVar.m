function [mu, varsigma] = ivmPosteriorMeanVar(model, X);

% IVMPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.

% IVM

D = size(model.y, 2);
numData = size(X, 1);

varsigma = zeros(numData, D);
mu = zeros(numData, D);

% Compute kernel for new point.
kX = kernCompute(model.kern, X, model.X(model.I, :))';

% COmpute diagonal of kernel for new point.
diagK = kernDiagCompute(model.kern, X);

if length(model.Sigma) > 1
  for i = 1:D
    if ~model.Sigma(i).robust
      Lk = model.Sigma(i).Linv*kX;
      Kinvk = model.Sigma(i).Linv'*Lk;
    else
      diagB = sparse(diag(sqrt(model.beta(model.I, i))));
      Lk = model.Sigma(i).Linv*diagB*kX;
      Kinvk = model.Sigma(i).Linv'*diagB*Lk;
    end
    for n = 1:numData
      varsigma(n, i) = diagK(n) - Lk(:, n)'*Lk(:, n); 
    end
  end
  mu(:, i) = Kinvk'*model.m(model.I, i);
else
  if ~model.Sigma.robust
    Lk = model.Sigma.Linv*kX;
    Kinvk = model.Sigma.Linv'*Lk;
  else
    diagB = sparse(diag(sqrt(model.beta(model.I, :))));
    Lk = model.Sigma.Linv*diagB*kX;
    Kinvk = model.Sigma.Linv'*diagB*Lk;
  end
  for n = 1:numData
    varsigma(n, :) = repmat(diagK(n) - Lk(:, n)'*Lk(:, n), 1, D);
  end
  mu = Kinvk'*model.m(model.I, :);
end
