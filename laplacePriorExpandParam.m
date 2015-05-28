function prior = laplacePriorExpandParam(prior, params)

% LAPLACEPRIOREXPANDPARAM Expand Laplace prior structure from param vector.

% GPMAT

prior.precision = params(1);
