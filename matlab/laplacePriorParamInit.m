function prior = laplacePriorParamInit(prior)

% LAPLACEPRIORPARAMINIT Laplace prior model's parameter initialisation.

% PRIOR

prior.precision = 1;

prior.transforms.index = [1];
prior.transforms.type = 'negLogLogit';
prior.nParams = 1;
