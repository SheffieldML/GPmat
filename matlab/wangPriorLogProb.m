function l = wangPriorLogProb(prior, x)

% WANGPRIORLOGPROB Log probability of Wang prior.

% SHEFFIELDML

% Compute log prior
D = length(x);
l = -prior.M*sum(log(x));
