function [params, names] = truncatedGammaPriorExtractParam(prior)

% TRUNCATEDGAMMAPRIOREXTRACTPARAM Extract params from truncated gamma prior structure.

% PRIOR

params = [prior.a prior.b];
if nargout > 1
  names = {'truncated gamma a', 'truncated gamma b'};
end
