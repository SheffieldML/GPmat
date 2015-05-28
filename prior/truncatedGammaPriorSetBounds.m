function prior = truncatedGammaPriorSetBounds(prior, bounds)

% TRUNCATEDGAMMAPRIORSETBOUNDS Set truncated gamma prior bounds.
%
% COPYRIGHT : Antti Honkela, 2013

% PRIOR

prior.lbound = bounds(1);
prior.ubound = bounds(2);
