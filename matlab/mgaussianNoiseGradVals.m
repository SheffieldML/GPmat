function [dlnZ_dmu, dlnZ_dvs] = mgaussianNoiseGradVals(noise, mu, varsigma, y)

% MGAUSSIANNOISEGRADVALS Gradient wrt mu and varsigma of log-likelihood for Variable variance Gaussian noise model.

% NOISE


D = size(y, 2);
nu = zeros(size(y));
dlnZ_dmu = zeros(size(y));
for i = 1:D
  nu(:, i) = 1./(noise.sigma2(i)+varsigma(:, i));
  dlnZ_dmu(:, i) = y(:, i) - mu(:, i) - noise.bias(i);
end
dlnZ_dmu = dlnZ_dmu.*nu;
dlnZ_dvs = -.5*nu+.5*dlnZ_dmu.*dlnZ_dmu;

% Remove missing values.
dlnZ_dmu(find(isnan(y))) = 0;
dlnZ_dvs(find(isnan(y))) = 0;
%/~
if any(isnan(dlnZ_dmu))
  warning('dlnZ_dmu is NaN')
end
if any(isnan(dlnZ_dvs))
  warning('dlnZ_dvs is NaN')
end
%~/