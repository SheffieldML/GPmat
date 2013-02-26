function prior = logisticNormalPriorExpandParam(prior, params)

% LOGISTICNORMALPRIOREXPANDPARAM Expand logistic-normal prior structure from params.
%
%	Description:
%	prior = logisticNormalPriorExpandParam(prior, params)
%

prior.mu = params(1);
prior.sd = params(2);
