function g = orderedGradX(X, Y, model, prior)

% ORDEREDGRADX Gradient wrt x of log-likelihood for Ordered categorical model.

% NOISE

% NOISE


if size(X, 1) > 1
  error('This function only takes one data-point');
end