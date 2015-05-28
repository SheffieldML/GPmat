function l = gaussianPriorLogProb(prior, x)

% GAUSSIANPRIORLOGPROB Log probability of Gaussian prior.

% PRIOR

% Compute log prior
l = -.5*sum(sum(prior.precision*x.*x + log(2*pi) - log(prior.precision)));
