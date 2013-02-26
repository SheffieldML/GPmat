function prior = normuniPriorExpandParam(prior, params)

% NORMUNIPRIOREXPANDPARAM Expand Normal uniform prior structure from param vector.

% SHEFFIELDML

prior.sigma = params(1);
prior.width = params(2);
