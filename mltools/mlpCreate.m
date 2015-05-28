function model = mlpCreate(inputDim, outputDim, options)

% MLPCREATE Multi-layer peceptron model.
% FORMAT
% DESC creates a structure for a multi-layer perceptron. For
% models with a single hidden layer it is a wrapper structure for NETLAB's
% multi-layer perceptron model.
% ARG inputDimension : dimension of input data.
% ARG outputDim : dimension of target data.
% ARG options : options structure. The structure contains the type
% of output 'activation function', the number of hidden units and
% the optimiser to be used. A set of default options are given by
% the file mlpOptions.
% RETURN model : model structure containing the neural network
% specified.
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2007
%
% SEEALSO : mlpOptions, mlp


% MLTOOLS

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
