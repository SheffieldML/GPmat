function [params, names] = laplacePriorExtractParam(prior)

% LAPLACEPRIOREXTRACTPARAM Extract params from Laplace prior structure.

% GPMAT

params = prior.precision;
if nargout > 1
  names = {'Laplace precision'};
end