function g = gaussianPriorGradient(prior, x)

% GAUSSIANPRIORGRADIENT Gradient wrt x of the log Gaussian prior.

% SHEFFIELDML

% Compute gradient of prior
g = -prior.precision*x;
