function g = noiseGradX(noise, mu, varsigma, dmu, dvs, y);

% NOISEGRADX Returns the gradient of the log-likelihood wrt x.

% NOISE

[dlnZ_dmu, dlnZ_dvs] = noiseGradVals(noise, mu, varsigma, y);
g = dlnZ_dmu*dmu' + dlnZ_dvs*dvs';
