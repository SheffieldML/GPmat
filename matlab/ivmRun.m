function model = ivmRun(XTrain, yTrain, kernelType, noiseType, ...
			selectionCriterion, dVal, options, kernelTieStructure, noiseTieStructure)

% IVMRUN Run ivm on a data set.

% IVM

% Intitalise IVM
model = ivm(XTrain, yTrain, kernelType, noiseType, selectionCriterion, dVal);
if nargin > 8 & ~isempty(noiseTieStructure);
  % Some of the noise parameters are constrained equal to each other
  model.noise = cmpndTieParameters(model.noise, noiseTieStructure);
end
if nargin > 7 & ~isempty(kernelTieStructure);
  % Some of the kernel parameters are constrained to equal each other.
  model.kern = cmpndTieParameters(model.kern, kernelTieStructure);
end

% Find the kernel parameters
model = ivmOptimise(model, options);

% Select data-points in an IVM with those kernel parameters.
model = ivmOptimiseIVM(model, options.display);
