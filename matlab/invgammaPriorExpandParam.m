function prior = invgammaPriorExpandParam(prior, params)

% INVGAMMAPRIOREXPANDPARAM Expand inverse gamma prior structure from params.

% IVM

prior.a = params(1);
prior.b = params(2);
