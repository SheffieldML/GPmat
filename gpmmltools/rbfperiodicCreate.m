function model = rbfperiodicCreate(inputDim, outputDim, options)

% RBFPERIODICCREATE Create a RBFPERIODIC model.
% This model is a periodic, single input, model. The model is based on
% constructing an RBF network in the two dimensional space given by x_1 =
% sin(theta) and x_2 = cos(theta).  This leads to basis functions of
% the form.
%
% phi(theta) = exp(-2/sigma2*sin^2(0.5*(theta - m)))
%
% SEEALSO : rbfperiodicKernCreate, rbfCreate
%
% FORMAT
% DESC creates a periodic radial basis function
% model structure given an options structure. 
% ARG inputDim : the input dimension of the model.
% ARG outputDim : the output dimension of the model.
% ARG options : an options structure that determines the form of the model.
% RETURN model : the model structure with the default parameters placed in.
%
% SEEALSO : rbfperiodicOptions, rbfperiodicParamInit, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2007

% MLTOOLS

if inputDim>1
  error(['You may only create a one dimensional input periodic RBF ' ...
         'model.'])
end
model.inputDim = inputDim;
model.outputDim = outputDim;

model.widthTransform.type = optimiDefaultConstraint('positive');
model.type = 'rbfperiodic';
model.hiddenDim = options.hiddenDim;

model.numParams = (inputDim+1)*options.hiddenDim + (options.hiddenDim + 1)*outputDim;
model = rbfperiodicParamInit(model);
