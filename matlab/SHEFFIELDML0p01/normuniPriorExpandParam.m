function prior = normuniPriorExpandParam(prior, params)

% NORMUNIPRIOREXPANDPARAM Expand Normal uniform prior structure from param vector.
%
%	Description:
%	prior = normuniPriorExpandParam(prior, params)
%

prior.sigma = params(1);
prior.width = params(2);
