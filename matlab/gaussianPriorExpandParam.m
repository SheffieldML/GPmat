function prior = gaussianPriorExpandParam(prior, params)

% GAUSSIANPRIOREXPANDPARAM Expand Gaussian prior structure from param vector.

% PRIOR

prior.precision = params(1);
