function l = mixturePriorLogProb(prior, x)

% MIXTUREPRIORLOGPROB Log probability of mixture prior.

% COPYRIGHT : Antti Honkela, 2013

% SHEFFIELDML

% Compute log prior

logprobs = cellfun(@(p) priorLogProb(p, x), prior.comp);
l = logsumexp(log(prior.weights) + logprobs);
