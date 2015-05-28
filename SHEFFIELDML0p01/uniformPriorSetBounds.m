function prior = uniformPriorSetBounds(prior, bounds)

% UNIFORMPRIORSETBOUNDS Set uniform prior bounds.
%
%	Description:
%	prior = uniformPriorSetBounds(prior, bounds)
%

prior.a = bounds(1);
prior.b = bounds(2);
