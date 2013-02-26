function prior = gammaPriorExpandParam(prior, params)

% GAMMAPRIOREXPANDPARAM Expand gamma prior structure from params.
%
%	Description:
%	prior = gammaPriorExpandParam(prior, params)
%

prior.a = params(1);
prior.b = params(2);
