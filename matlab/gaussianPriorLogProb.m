function l = gaussianPriorLogProb(prior, x)

% GAUSSIANPRIORLOGPROB Log probability of Gaussian prior.

% IVM

% Compute log prior
l = -.5*(prior.precision*x*x' + log(2*pi) - log(prior.precision));
