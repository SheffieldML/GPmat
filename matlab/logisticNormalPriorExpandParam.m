function prior = logisticNormalPriorExpandParam(prior, params)

% LOGISTICNORMALPRIOREXPANDPARAM Expand logistic-normal prior structure from params.

% SHEFFIELDML

prior.mu = params(1);
prior.sd = params(2);
