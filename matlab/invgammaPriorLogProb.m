function l = invgammaPriorLogProb(prior, x)

% INVGAMMAPRIORLOGPROB Log probability of inverse gamma prior.

% PRIOR

% PRIOR


% Compute log prior
D = length(x);
l = D*prior.a*log(prior.b)-D*gammaln(prior.a)-(prior.a+1)*sum(log(x))-prior.b*sum(1./x);
