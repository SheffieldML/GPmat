function [dlnZ_dmu, dlnZ_dvs] = ngaussNoiseGradVals(noise, mu, varsigma, y)

% NGAUSSNOISEGRADVALS Gradient wrt mu and varsigma of log-likelihood for noiseless Gaussian noise model.

% IVM

[dlnZ_dmu, dlnZ_dvs] = gaussianNoiseGradVals(noise, mu, varsigma, y);