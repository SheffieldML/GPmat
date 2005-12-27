function [mu, varsigma] = fgplvmDynamicsPosteriorMeanVar(model, X);

% FGPLVMDYNAMICSPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.

% FGPLVM


varsigma = zeros(size(X, 1), model.q);
mu = zeros(size(X, 1), model.q);

% Compute kernel for new point.
kX = kernCompute(model.dynamics.kern, X, model.X(1:end-1, :))';
  
invKk = model.dynamics.invK_uu*kX;

if nargout > 1
  % Compute diagonal of kernel for new point.
  diagK = kernDiagCompute(model.dynamics.kern, X);
  for n = 1:size(X, 1)
    varsigma(n, :) = repmat(diagK(n) - kX(:, n)'*invKk(:, n), 1, model.q);
    if varsigma(n, 1) < 0
      warning('Varsigma less than zero');
    end
  end
end
mu = invKk'*model.X(2:end, :);
