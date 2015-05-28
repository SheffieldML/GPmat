function [params, names] = uniformPriorExtractParam(prior)

% UNIFORMPRIOREXTRACTPARAM Extract params from uniform prior structure.

% GPMAT

params = [];
if nargout > 1
  names = {};
end
