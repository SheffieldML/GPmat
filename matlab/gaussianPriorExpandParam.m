function prior = gaussianPriorExpandParam(prior, params)

% GAUSSIANPRIOREXPANDPARAM Expand Gaussian prior structure from param vector.

% IVM

prior.precision = params(1);
