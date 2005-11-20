function prior = wangPriorParamInit(prior)

% WANGPRIORPARAMINIT Wang prior model's parameter initialisation.

% PRIOR

prior.M = 1;

prior.transforms.index = [1];
prior.transforms.type = 'negLogLogit';
prior.nParams = 1;
