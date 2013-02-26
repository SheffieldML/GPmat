function l = wangPriorLogProb(prior, x)

% WANGPRIORLOGPROB Log probability of Wang prior.
%
%	Description:
%	l = wangPriorLogProb(prior, x)
%

% Compute log prior
D = length(x);
l = -prior.M*sum(log(x));
