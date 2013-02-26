function prior = logisticNormalPriorSetBounds(prior, bounds)

% LOGISTICNORMALPRIORSETBOUNDS Set logistic-normal prior bounds.

% SHEFFIELDML

prior.a = bounds(1);
prior.b = bounds(2);
