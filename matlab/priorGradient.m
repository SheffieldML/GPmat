function g = priorGradient(prior, params)

% PRIORGRADIENTPARAM Gradient of the prior with respect to its variables

% PRIOR

g = feval([prior.type 'PriorGradient'], prior, params);