function prior = gammaPriorExpandParam(prior, params)

% GAMMAPRIOREXPANDPARAM Expand gamma prior structure from params.

% IVM

prior.a = params(1);
prior.b = params(2);
