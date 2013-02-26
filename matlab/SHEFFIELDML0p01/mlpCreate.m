function model = mlpCreate(inputDim, outputDim, options)

% MLPCREATE Multi-layer peceptron model.
%
%	Description:
%
%	MODEL = MLPCREATE(INPUTDIMENSION, OUTPUTDIM, OPTIONS) creates a
%	structure for a multi-layer perceptron. For models with a single
%	hidden layer it is a wrapper structure for NETLAB's multi-layer
%	perceptron model.
%	 Returns:
%	  MODEL - model structure containing the neural network specified.
%	 Arguments:
%	  INPUTDIMENSION - dimension of input data.
%	  OUTPUTDIM - dimension of target data.
%	  OPTIONS - options structure. The structure contains the type of
%	   output 'activation function', the number of hidden units and the
%	   optimiser to be used. A set of default options are given by the
%	   file mlpOptions.
%	
%
%	See also
%	MLPOPTIONS, MLP


%	Copyright (c) 2005, 2006, 2007 Neil D. Lawrence



if length(options.hiddenDim) == 1
  % Can use NETLAB implementation.
  model = mlp(inputDim, options.hiddenDim, outputDim, ...
              options.activeFunc);
  model.numParams = model.nwts;
  model.hiddenDim = options.hiddenDim;
  model.inputDim = inputDim;
  model.outputDim = outputDim;
else
  error('Multiple hidden layer mlp not yet implemented')
  model.type = 'mlp';
  model.hiddenDim = options.hiddenDim;
  model.inputDim = inputDim;
  model.outputDim = outputDim;
  model.numParams = (model.inputDim+1)*model.hiddenDim(1);
  for i = 1:length(model.hiddenDim)-1
    model.numParams = model.numParams + (model.hiddenDim(i)+1)*model.hiddenDim(i+1);
  end
  model.numParams = model.numParams + (model.hiddenDim(end)+1)*model.outputDim;

  activationFunctions = {'linear', 'logistic', 'softmax'};

  if sum(strcmp(options.activeFunc, activationFunctions)) == 0
    error('Undefined output function.');
  else
    model.outfn = options.activeFunc;
  end
  model = mlpParamInit(model);
end

model.optimiser = options.optimiser;
