function [params, names] = uniformPriorExtractParam(prior)

% UNIFORMPRIOREXTRACTPARAM Extract params from uniform prior structure.

% SHEFFIELDML

params = [];
if nargout > 1
  names = {};
end
