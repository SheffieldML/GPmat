function l = gammaPriorLogProb(prior, x)

% GAMMAPRIORLOGPROB Log probability of Gamma prior.

% PRIOR

% Compute log prior
D = length(x);
l = D*prior.a*log(prior.b)-D*gammaln(prior.a)+(prior.a-1)*sum(log(x))-prior.b*sum(x);
