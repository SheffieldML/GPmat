function [params, names] = gaussianPriorExtractParam(prior)

% GAUSSIANPRIOREXTRACTPARAM Extract params from Gaussian prior structure.

% PRIOR

params = prior.precision;
if nargout > 1
  names = {'Gaussian precision'};
end
