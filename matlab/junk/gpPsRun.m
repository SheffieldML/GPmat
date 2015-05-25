function models = gpPsRun(X, y, kernelType, noiseType, theta, prior, display, iters)

% GPPSRUN Run a GP on point-set data.

% PSIVM

models = psivm(X, y, kernelType, noiseType, 'none');

models.lntheta = log(thetaConstrain(theta));
models = psivmInit(models);

models = psivmOptimiseKernel(models, prior, display, iters);
