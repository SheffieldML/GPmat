function l = priorLogProb(prior, x)

% PRIORLOGPROB Log probability of Gaussian prior.

% PRIOR

% PRIOR


% Compute log prior
l = feval([prior.type 'PriorLogProb'], prior, x);
