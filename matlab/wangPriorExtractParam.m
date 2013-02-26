function [params, names] = wangPriorExtractParam(prior)

% WANGPRIOREXTRACTPARAM Extract params from Wang prior structure.

% SHEFFIELDML

params = [prior.M];
if nargout > 1
  names = {'M'};
end