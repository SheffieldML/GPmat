function g = laplacePriorGradient(prior, x)

% LAPLACEPRIORGRADIENT Gradient wrt x of the log Laplace prior.

% PRIOR

% Compute gradient of prior
g = -prior.precision*sign(x);
