function model = ivmRun(XTrain, yTrain, kernelType, noiseType, ...
			selectionCriterion, dVal, prior, display, innerIters, ...
			   outerIters, kernelTieStructure, noiseTieStructure)

% IVMRUN Run ivm on a data set.

% IVM

% Intitalise IVM
model = ivm(XTrain, yTrain, kernelType, noiseType, selectionCriterion, dVal);
if nargin > 11 & ~isempty(noiseTieStructure);
  % Some of the noise parameters are constrained equal to each other
  model.noise = cmpndTieParameters(model.noise, noiseTieStructure);
end
if nargin > 10 & ~isempty(kernelTieStructure);
  % Some of the kernel parameters are constrained to equal each other.
  model.kern = cmpndTieParameters(model.kern, kernelTieStructure);
end

% Find the kernel parameters
model = ivmOptimise(model, prior, display, innerIters, outerIters);

% Select data-points in an IVM with those kernel parameters.
model = ivmOptimiseIVM(model, display);
