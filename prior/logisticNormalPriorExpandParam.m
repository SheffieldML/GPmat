function prior = logisticNormalPriorExpandParam(prior, params)

% LOGISTICNORMALPRIOREXPANDPARAM Expand logistic-normal prior structure from params.

% PRIOR

prior.mu = params(1);
prior.sd = params(2);
