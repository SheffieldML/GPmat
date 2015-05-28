function [params, names] = gammaPriorExtractParam(prior)

% GAMMAPRIOREXTRACTPARAM Extract params from gamma prior structure.

% GPMAT

params = [prior.a prior.b];
if nargout > 1
  names = {'gamma a', 'gamma b'};
end