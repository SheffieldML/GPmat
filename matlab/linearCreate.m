function model = linearCreate(inputDim, outputDim, options)

% LINEARCREATE Create a linear model.

% MLTOOLS

model.type = 'linear';
model.activeFunc = options.activeFunc; 
model.inputDim = inputDim;
model.outputDim = outputDim;
model.numParams = (inputDim + 1)*outputDim;


model.W = randn(inputDim, outputDim)/sqrt(inputDim+1);
model.b = randn(1, outputDim)/sqrt(inputDim+1);
