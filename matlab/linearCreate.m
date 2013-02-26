function model = linearCreate(inputDim, outputDim, options)

% LINEARCREATE Create a linear model.

% SHEFFIELDML

model.type = 'linear';
model.activeFunc = options.activeFunc; 
model.inputDim = inputDim;
model.outputDim = outputDim;
model.numParams = (inputDim + 1)*outputDim;

model = linearParamInit(model);
