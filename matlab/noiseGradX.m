function g = noiseGradX(noise, mu, varsigma, dmu, dvs, y);

% NOISEGRADX Returns the gradient of the log-likelihood wrt x.

[dlnZ_dmu, dlnZ_dvs] = feval([noise.type 'NoiseGradVals'], noise, mu, varsigma, y);
g = dlnZ_dmu*dmu' + dlnZ_dvs*dvs';
