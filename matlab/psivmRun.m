function models = psivmRun(XTrain, yTrain, kernelType, noiseType, ...
			selectionCriterion, dVal, prior, display, innerIters, ...
			   outerIters)

% PSIVMRUN Run point set ivm on a data set.

% PSIVM

% Intitalise IVM
models = psivm(XTrain, yTrain, kernelType, noiseType, selectionCriterion, dVal);
models = psivmOptimise(models, prior, display, innerIters, outerIters);
