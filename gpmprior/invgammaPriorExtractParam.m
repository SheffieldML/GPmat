function [params, names] = invgammaPriorExtractParam(prior)

% INVGAMMAPRIOREXTRACTPARAM Extract params from inverse gamma prior structure.

% PRIOR

% PRIOR


params = [prior.a prior.b];
if nargout > 1
  names = {'inv gamma a', 'inv gamma b'};
end
