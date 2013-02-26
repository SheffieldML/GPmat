function prior = gammaPriorExpandParam(prior, params)

% GAMMAPRIOREXPANDPARAM Expand gamma prior structure from params.

% SHEFFIELDML

prior.a = params(1);
prior.b = params(2);
