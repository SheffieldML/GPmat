function prior = normuniPriorParamInit(prior)

% NORMUNIPRIORPARAMINIT Normal uniform prior model's parameter initialisation.

% PRIOR

prior.sigma = 0.1;
prior.width = 2;
prior.transforms.index = [1 2];
prior.transforms.type = 'negLogLogit';
prior.nParams = 2;
