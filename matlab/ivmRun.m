function model = ivmRun(XTrain, yTrain, kernelType, noiseType, ...
			selectionCriterion, dVal, prior, display, innerIters, ...
			   outerIters, tieStructure)

% IVMRUN Run ivm on a data set.

% IVM

% Intitalise IVM
model = ivm(XTrain, yTrain, kernelType, noiseType, selectionCriterion, dVal);

if nargin > 10
  % Some of the kernel parameters are constrained to equal each other.
  model.kern = cmpndTieParameters(model.kern, tieStructure);
end

% Find the kernel parameters
model = ivmOptimise(model, prior, display, innerIters, outerIters);

% Select data-points in an IVM with those kernel parameters.
model = ivmOptimiseIVM(model, display);
