function model = kbrCreate(inputDim, outputDim, options)

% KBRCREATE Create a KBR model.
% The kernel based regression model is simply a model for least
% squares regression in a kernel feature space. Any kernel from the KERN
% toolbox can be specified. The model was developed for providing kernel
% based back constraints in the GP-LVM. Please consider using a Gaussian
% process model (through the GP toolbox) if you are interested in the
% model for regression.
%
% FORMAT
% DESC creates a kernel based regression
%  model structure given an options structure. 
% ARG options : an options structure that determines the form of the model.
% RETURN model : the model structure with the default parameters placed in.
%
% SEEALSO : kbrOptions, kbrParamInit, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% MLTOOLS


model.type = 'kbr';
model.inputDim = inputDim;
model.outputDim = outputDim;
model.numData = size(options.X, 1);
model.numParams = (model.numData + 1)*outputDim;
model.X = options.X;
if isstruct(options.kern)
  model.kern = options.kern;
else
  model.kern = kernCreate(options.X, options.kern);
end

model.K = kernCompute(model.kern, options.X);
model = kbrParamInit(model);
