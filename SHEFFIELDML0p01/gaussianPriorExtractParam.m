function [params, names] = gaussianPriorExtractParam(prior)

% GAUSSIANPRIOREXTRACTPARAM Extract params from Gaussian prior structure.
%
%	Description:
%	[params, names] = gaussianPriorExtractParam(prior)
%

params = prior.precision;
if nargout > 1
  names = {'Gaussian precision'};
end