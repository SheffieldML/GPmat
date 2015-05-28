function g = normuniPriorGradient(prior, x)

% NORMUNIPRIORGRADIENT Gradient wrt x of the log normal uniform prior.

% PRIOR

% Compute gradient of prior
u = (x+prior.width/2)/prior.sigma;
uprime = u - prior.width/prior.sigma;

B1 = gaussOverDiffCumGaussian(u, uprime, 1);
B2 = gaussOverDiffCumGaussian(u, uprime, 2);
g = (B1-B2)/prior.sigma;

