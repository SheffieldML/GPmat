function model = mlpCreate(inputDim, outputDim, options)

% MLPCREATE Wrapper for NETLAB's mlp `net'.
% FORMAT
% DESC creates a structure used as a wrapper structure for NETLAB's
% multi-layer perceptron model.
% ARG inputDimension : dimension of input data.
% ARG outputDim : dimension of target data.
% ARG options : options structure. The structure contains the type
% of output 'activation function', the number of hidden units and
% the optimiser to be used. A set of default options are given by
% the file mlpOptions.
% RETURN model : model structure containing the neural networks
% specified.
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% SEEALSO : mlpOptions, mlp


% MLTOOLS

model = mlp(inputDim, options.hiddenDim, outputDim, options.activeFunc);
model.optimiser = options.optimiser;
model.numParams = model.nwts;
model.inputDim = inputDim;
model.outputDim = outputDim;
