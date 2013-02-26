function [params, names] = wangPriorExtractParam(prior)

% WANGPRIOREXTRACTPARAM Extract params from Wang prior structure.
%
%	Description:
%	[params, names] = wangPriorExtractParam(prior)
%

params = [prior.M];
if nargout > 1
  names = {'M'};
end