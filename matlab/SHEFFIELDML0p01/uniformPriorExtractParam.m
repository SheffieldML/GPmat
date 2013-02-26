function [params, names] = uniformPriorExtractParam(prior)

% UNIFORMPRIOREXTRACTPARAM Extract params from uniform prior structure.
%
%	Description:
%	[params, names] = uniformPriorExtractParam(prior)
%

params = [];
if nargout > 1
  names = {};
end
