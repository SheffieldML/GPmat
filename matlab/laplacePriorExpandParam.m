function prior = laplacePriorExpandParam(prior, params)

% LAPLACEPRIOREXPANDPARAM Expand Laplace prior structure from param vector.

% SHEFFIELDML

prior.precision = params(1);
