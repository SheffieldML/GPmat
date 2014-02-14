function prior = uniformPriorSetBounds(prior, bounds)

% UNIFORMPRIORSETBOUNDS Set uniform prior bounds.

% GPMAT

prior.a = bounds(1);
prior.b = bounds(2);
