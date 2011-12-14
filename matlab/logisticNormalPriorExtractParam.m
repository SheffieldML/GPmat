function [params, names] = logisticNormalPriorExtractParam(prior)

% LOGISTICNORMALPRIOREXTRACTPARAM Extract params from logistic-normal prior structure.

% PRIOR

params = [prior.mu, prior.sd, prior.a, prior.b];
if nargout > 1
  names = {'logistic-normal mean', 'logistic-normal std', ...
	   'logistic-normal a', 'logistic-normal b'};
end
