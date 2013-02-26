function prior = logisticNormalPriorSetBounds(prior, bounds)

% LOGISTICNORMALPRIORSETBOUNDS Set logistic-normal prior bounds.
%
%	Description:
%	prior = logisticNormalPriorSetBounds(prior, bounds)
%

prior.a = bounds(1);
prior.b = bounds(2);
