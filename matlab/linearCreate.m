function model = linearCreate(inputDim, outputDim)

% LINEARCREATE Create a linear model.

model.type = 'linear';
model.inputDim = inputDim;
model.outputDim = outputDim;
model.numParams = (inputDim + 1)*outputDim;


model.W = randn(inputDim, outputDim)/sqrt(inputDim+1);
model.b = randn(1, outputDim)/sqrt(inputDim+1);
