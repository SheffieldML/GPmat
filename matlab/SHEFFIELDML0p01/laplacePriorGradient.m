function g = laplacePriorGradient(prior, x)

% LAPLACEPRIORGRADIENT Gradient wrt x of the log Laplace prior.
%
%	Description:
%	g = laplacePriorGradient(prior, x)
%

% Compute gradient of prior
g = -prior.precision*sign(x);
