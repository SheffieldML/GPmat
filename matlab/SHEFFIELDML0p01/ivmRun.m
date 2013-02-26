function model = ivmRun(XTrain, yTrain, options, kernelTieStructure, noiseTieStructure)

% IVMRUN Run the IVM on a given data set.
%
%	Description:
%
%	IVMRUN(X, Y, OPTIONS) takes an input data matrix, a target data
%	matrix and an options structure (as defined by ivmOptions). It then
%	creates and optimises an IVM using ivmOptimise.
%	 Arguments:
%	  X - the input training data.
%	  Y - the target training data.
%	  OPTIONS - options structure (as defined by ivmOptions).
%
%	IVMRUN(X, Y, OPTIONS, KERNELTIE, NOISETIE) is as above but includes
%	the facility to 'tie' some of the kernel an noise parameters
%	together.
%	 Arguments:
%	  X - the input training data.
%	  Y - the target training data.
%	  OPTIONS - options structure (as defined by ivmOptions).
%	  KERNELTIE - structure which speficies which kernel parameters are
%	   forced to be equal to each other (see cmpndTieParameters).
%	  NOISETIE - structure which speficies which noise parameters are
%	   forced to be equal to each other (see cmpndTieParameters).
%	
%
%	See also
%	IVMCREATE, IVMOPTIMISE, IVMOPTIONS, CMPNDTIEPARAMETERS


%	Copyright (c) 2004, 2005, 2006, 2007 Neil D. Lawrence


% Intitalise IVM
model = ivmCreate(size(XTrain, 1), size(yTrain, 2), XTrain, yTrain, options);

if nargin > 8 & ~isempty(noiseTieStructure);
  % Some of the noise parameters are constrained equal to each other
  model.noise = cmpndTieParameters(model.noise, noiseTieStructure);
end
if nargin > 7 & ~isempty(kernelTieStructure);
  % Some of the kernel parameters are constrained to equal each other.
  model.kern = cmpndTieParameters(model.kern, kernelTieStructure);
end

% Optimise the parameters
model = ivmOptimise(model, options);

% Select data-points in an IVM with the given parameters.
model = ivmOptimiseIvm(model, options.display);
