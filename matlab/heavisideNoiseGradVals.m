function [dlnZ_dmu, dlnZ_dvs] = heavisideNoiseGradVals(noise, mu, varsigma, y)

% HEAVISIDENOISEGRADVALS Gradient wrt mu and varsigma of log-likelihood for heaviside noise model.

% IVM

D = size(mu, 2);
c = y./sqrt(varsigma);
u = zeros(size(c));
dlnZ_dmu = zeros(size(c));
for i = 1:D
  u(:, i) = c(:, i).*(mu(:, i) + noise.bias(i));
  dlnZ_dmu(:, i) = c(:, i).*ngaussian(u(:, i))...
            ./(cumGaussian(u(:, i))+ ...
               noise.eta/(1-2*noise.eta));
end
dlnZ_dvs = -.5*c.*u.*dlnZ_dmu;
