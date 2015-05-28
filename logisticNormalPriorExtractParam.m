function [params, names] = logisticNormalPriorExtractParam(prior)

% LOGISTICNORMALPRIOREXTRACTPARAM Extract params from logistic-normal prior structure.

% PRIOR

params = [prior.mu, prior.sd];
if nargout > 1
  names = {'logistic-normal mean', 'logistic-normal std'};
end
