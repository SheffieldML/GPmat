function g = gnmGradX(X, Y, model, prior)

% GNMGRADX Gradient wrt x of log-likelihood for gap noise model.

% IVM

if size(X, 1) > 1
  error('This function only takes one data-point');
end