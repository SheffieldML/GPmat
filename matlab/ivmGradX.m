function g = ivmGradX(model, x, y);

% IVMGRADX Returns the gradient of the log-likelihood wrt x.

[mu, varsigma] = ivmPosteriorMeanVar(x, model);
[dmu, dvs] = ivmPosteriorGradMeanVar(x, model);

g = noiseGradX(model.noise, mu, varsigma, dmu, dvs, y);