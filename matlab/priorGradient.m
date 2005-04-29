function g = priorGradient(prior, params)

% PRIORGRADIENT Gradient of the prior with respect to its variables

% PRIOR
fhandle = str2func([prior.type 'PriorGradient']);
g = fhandle(prior, params);