function prior = priorParamInit(prior)

% PRIORPARAMINIT Prior model's parameter initialisation.

% IVM

prior = feval([prior.type 'PriorParamInit'], prior);
