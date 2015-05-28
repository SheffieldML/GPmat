function prior = laplacePriorExpandParam(prior, params)

% LAPLACEPRIOREXPANDPARAM Expand Laplace prior structure from param vector.

% PRIOR

prior.precision = params(1);
