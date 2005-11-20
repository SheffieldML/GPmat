function prior = priorCreate(type)

% PRIORCREATE Create a prior structure given a type.

% PRIOR

prior.type = type;
prior = priorParamInit(prior);
