function [nu, dlnZ_dmu, dlnZ_dvs] = cmpndGradVals(noise, mu, varsigma, y)

% CMPNDGRADX Gradient wrt x of log-likelihood for compound noise model.

% IVM

if size(X, 1) > 1
  error('This function only takes one data-point');
end