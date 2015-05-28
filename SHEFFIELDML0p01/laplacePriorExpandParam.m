function prior = laplacePriorExpandParam(prior, params)

% LAPLACEPRIOREXPANDPARAM Expand Laplace prior structure from param vector.
%
%	Description:
%	prior = laplacePriorExpandParam(prior, params)
%

prior.precision = params(1);
