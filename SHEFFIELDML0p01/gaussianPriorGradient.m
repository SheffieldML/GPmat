function g = gaussianPriorGradient(prior, x)

% GAUSSIANPRIORGRADIENT Gradient wrt x of the log Gaussian prior.
%
%	Description:
%	g = gaussianPriorGradient(prior, x)
%

% Compute gradient of prior
g = -prior.precision*x;
