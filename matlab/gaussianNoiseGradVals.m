function [dlnZ_dmu, dlnZ_dvs] = gaussianNoiseGradVals(noise, mu, varsigma, y)

% GAUSSIANNOISEGRADVALS Gradient wrt mu and varsigma of log-likelihood for gaussian noise model.

% IVM

nu = 1./(varsigma+noise.sigma2);
dlnZ_dmu = (y-(mu+noise.bias)).*nu;
dlnZ_dvs = -.5*nu+.5*dlnZ_dmu.*dlnZ_dmu;
