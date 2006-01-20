function model = kbrCreate(inputDim, outputDim, options)

% KBRCREATE Create a kernel based regression model.

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
model.A = randn(model.numData, outputDim)/sqrt(model.numData+1);
model.bias = randn(1, outputDim)/sqrt(model.numData+1);
