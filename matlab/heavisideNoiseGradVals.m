function [dlnZ_dmu, dlnZ_dvs] = heavisideNoiseGradVals(noise, mu, varsigma, y)

% HEAVISIDENOISEGRADVALS Gradient wrt mu and varsigma of log-likelihood for heaviside noise model.

% IVM

c = y./sqrt(varsigma);
u = zeros(size(c));
dlnZ_dmu = zeros(size(c));
u = c.*(mu+noise.bias);
dlnZ_dmu = c.*ngaussian(u)./(cummGaussian(u) + noise.eta./(1-2*noise.eta)); 

dlnZ_dvs = -.5*c.*u.*dlnZ_dmu;
