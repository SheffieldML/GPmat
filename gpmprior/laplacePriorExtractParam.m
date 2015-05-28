function [params, names] = laplacePriorExtractParam(prior)

% LAPLACEPRIOREXTRACTPARAM Extract params from Laplace prior structure.

% PRIOR

params = prior.precision;
if nargout > 1
  names = {'Laplace precision'};
end
