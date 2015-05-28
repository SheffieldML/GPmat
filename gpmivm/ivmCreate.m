function model = ivmCreate(inputDim, outputDim, X, y, options);

% IVMCREATE Create a IVM model with the IVM sparse approximaiton.
% The IVM stands for Informative Vector Machine. The IVM is a
% sparse Gaussian process approximation which uses information
% theoretic criteria for selection of an active set. Unlike most
% sparse Gaussian process methods (but like the SVM) the IVM is a
% compression scheme. In other words the IVM can be recreated using
% only those data points considered to be 'informative vectors'. We
% refer to these points as the active set.
%
% FORMAT
% DESC creates a Gaussian process model structure using the IVM
% approximation. The default parameter settings as specified by the
% options vector. This function replaces the deprecated ivm function.
% ARG q : input data dimension.
% ARG d : the number of processes (i.e. output data dimension).
% ARG X : the input data matrix.
% ARG y : the target (output) data.
% ARG options : options structure as defined by ivmOptions.m.
% RETURN model : model structure containing the IVM.
%
% SEEALSO : ivmOptions, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2007

% IVM


model.type = 'ivm';
model.terminate = 0;
model.epUpdate = 0;

model.d = options.numActive;

model.X = X;
model.y = y;

model.m = [];
model.beta = [];

model.nu = zeros(size(y));
model.g = zeros(size(y));

if isstruct(options.kern) 
  model.kern = options.kern;
else
  model.kern = kernCreate(model.X, options.kern);
end

model.varSigma = zeros(size(y));
model.mu = zeros(size(y));

model.I = [];
model.J = [];

if isstruct(options.noise) 
  model.noise = options.noise;
else
  model.noise = noiseCreate(options.noise, model.y);
end

if model.noise.spherical
  model.Sigma.M = [];
  model.Sigma.L = [];
else
  for i = 1:size(y, 2)
    model.Sigma(i).M = [];
    model.Sigma(i).L = [];
  end
end
model.selectionCriterion = options.selectionCriterion;


switch options.selectionCriterion
 case 'none'
  numData = size(X, 1);
  model.I = (1:numData);
 otherwise
  
end
