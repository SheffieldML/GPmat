function prior = priorParamInit(prior)

% PRIORPARAMINIT Prior model's parameter initialisation.

% PRIOR

% PRIOR


prior = feval([prior.type 'PriorParamInit'], prior);
