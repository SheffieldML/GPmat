function g = priorGradient(prior, params)

% PRIORGRADIENTPARAM Gradient of the prior with respect to its variables

g = feval([prior.type 'PriorGradient'], prior, params);