function l = laplacePriorLogProb(prior, x)

% LAPLACEPRIORLOGPROB Log probability of Laplace prior.

% SHEFFIELDML

% Compute log prior
l = -prior.precision*sum(abs(x)) - log(2) + log(prior.precision);
