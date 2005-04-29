function l = priorLogProb(prior, x)

% PRIORLOGPROB Log probability of Gaussian prior.

% PRIOR

% Compute log prior
fhandle = str2func([prior.type 'PriorLogProb']);
l = fhandle(prior, x);
