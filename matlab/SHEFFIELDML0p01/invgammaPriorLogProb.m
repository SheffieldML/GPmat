function l = invgammaPriorLogProb(prior, x)

% INVGAMMAPRIORLOGPROB Log probability of inverse gamma prior.
%
%	Description:
%	l = invgammaPriorLogProb(prior, x)
%



% Compute log prior
D = length(x);
l = D*prior.a*log(prior.b)-D*gammaln(prior.a)-(prior.a+1)*sum(log(x))-prior.b*sum(1./x);
