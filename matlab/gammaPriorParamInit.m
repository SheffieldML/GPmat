function prior = gammaPriorParamInit(prior)

% GAMMAPRIORPARAMINIT Gamma prior model's parameter initialisation.

% IVM

prior.a = 1e-6;
prior.b = 1e-6;

prior.transforms.index = [1 2];
prior.transforms.type = 'negLogLogit';
prior.nParams = 2;
