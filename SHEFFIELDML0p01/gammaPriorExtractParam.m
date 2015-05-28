function [params, names] = gammaPriorExtractParam(prior)

% GAMMAPRIOREXTRACTPARAM Extract params from gamma prior structure.
%
%	Description:
%	[params, names] = gammaPriorExtractParam(prior)
%

params = [prior.a prior.b];
if nargout > 1
  names = {'gamma a', 'gamma b'};
end