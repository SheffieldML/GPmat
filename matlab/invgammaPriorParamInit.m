function prior = invgammaPriorParamInit(prior)

% INVGAMMAPRIORPARAMINIT Inverse gamma prior model's parameter initialisation.

% PRIOR


prior.a = 1e-6;
prior.b = 1e-6;

prior.transforms.index = [1 2];
prior.transforms.type = 'negLogLogit';
prior.nParams = 2;
