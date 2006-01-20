function model = mlpCreate(inputDim, outputDim, options)

% MLPCREATE Wrapper for NETLAB's mlp `net'.

% MLTOOLS

model = mlp(inputDim, options.hiddenDim, outputDim, options.activeFunc);
model.numParams = model.nwts;
model.inputDim = inputDim;
model.outputDim = outputDim;
