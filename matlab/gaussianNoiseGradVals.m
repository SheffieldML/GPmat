function [dlnZ_dmu, dlnZ_dvs] = gaussianNoiseGradVals(noise, mu, varsigma, y)

% GAUSSIANNOISEGRADVALS Gradient wrt mu and varsigma of log-likelihood for gaussian noise model.

% IVM


D = size(y, 2);
nu = 1./(noise.sigma2+varsigma);
dlnZ_dmu = zeros(size(nu));
for i = 1:D
  dlnZ_dmu(:, i) = y(:, i) - mu(:, i) - noise.bias(i);
end
dlnZ_dmu = dlnZ_dmu.*nu;
dlnZ_dvs = -.5*nu+.5*dlnZ_dmu.*dlnZ_dmu;
