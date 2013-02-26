function [params, names] = laplacePriorExtractParam(prior)

% LAPLACEPRIOREXTRACTPARAM Extract params from Laplace prior structure.
%
%	Description:
%	[params, names] = laplacePriorExtractParam(prior)
%

params = prior.precision;
if nargout > 1
  names = {'Laplace precision'};
end