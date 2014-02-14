function l = gaussian2PriorLogProb(prior, x)

% GAUSSIAN2PRIORLOGPROB Log probability of Gaussian prior.
% COPYRIGHT: Neil D. Lawrence 2004
% MODIFICATIONS: Andreas C. Damianou, 2013
% SHEFFIELDML

% Compute log prior
x_c = x - prior.mean;
l = length(prior.mean) * (-.5*sum(sum(prior.precision*x_c.*x_c + log(2*pi) - log(prior.precision))));
