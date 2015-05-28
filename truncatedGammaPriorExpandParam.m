function prior = truncatedGammaPriorExpandParam(prior, params)

% TRUNCATEDGAMMAPRIOREXPANDPARAM Expand truncated gamma prior structure from params.

% PRIOR

prior.a = params(1);
prior.b = params(2);
