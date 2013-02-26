function l = uniformPriorLogProb(prior, x)

% UNIFORMPRIORLOGPROB Log probability of uniform prior.
%
%	Description:
%	l = uniformPriorLogProb(prior, x)
%

% Compute log prior

l = sum(log((x >= prior.a) .* (x <= prior.b) * 1 / (prior.b - prior.a)));
