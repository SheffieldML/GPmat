function g = gaussian2PriorGradient(prior, x)

% GAUSSIAN2PRIORGRADIENT Gradient wrt x of the log Gaussian prior.
% COPYRIGHT: Neil D. Lawrence 2004
% MODIFICATIONS: Andreas C. Damianou, 2013
% SHEFFIELDML

% Compute gradient of prior
g = -(prior.precision*length(x)).*(x - prior.mean);
