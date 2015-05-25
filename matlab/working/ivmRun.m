function model = ivmRun(XTrain, yTrain, kernelType, noiseType, ...
			selectionCriterion, dVal, prior, display, innerIters, ...
			   outerIters)

% IVMRUN Run ivm on a data set.

% IVM

% Intitalise IVM
model = ivm(XTrain, yTrain, kernelType, noiseType, selectionCriterion, dVal);

% Find the kernel parameters
model = ivmOptimise(model, prior, display, innerIters, outerIters);

% Do final point inclusion
model = ivmOptimiseIVM(model, display);
