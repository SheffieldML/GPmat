function g = ivmGradX(model, x, y);

% IVMGRADX Returns the gradient of the log-likelihood wrt x.

[mu, varsigma] = ivmPosteriorMeanVar(model, x);
[dmu, dvs] = ivmPosteriorGradMeanVar(model, x);

g = noiseGradX(model.noise, mu, varsigma, dmu, dvs, y);