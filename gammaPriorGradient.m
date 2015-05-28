function g = gammaPriorGradient(prior, x)

% GAMMAPRIORGRADIENT Gradient wrt x of the gamma prior.

% GPMAT

% Compute gradient of prior
g = (prior.a-1)./x-prior.b;
