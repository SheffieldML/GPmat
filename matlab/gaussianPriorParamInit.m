function prior = gaussianPriorParamInit(prior)

% GAUSSIANPRIORPARAMINIT Gaussian prior model's parameter initialisation.

% IVM

prior.precision = 1;

prior.transforms.index = [1];
prior.transforms.type = 'negLogLogit';
prior.nParams = 1;
