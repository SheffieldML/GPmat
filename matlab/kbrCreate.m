function model = kbrCreate(inputDim, outputDim, kern, X)

% KBRCREATE Create a kernel based regression model.

model.type = 'kbr';
model.inputDim = inputDim;
model.outputDim = outputDim;
model.numData = size(X, 1);
model.numParams = (model.numData + 1)*outputDim;
model.X = X;
if isstruct(kern)
  model.kern = kern;
else
  model.kern = kernCreate(X, kern);
end

model.K = kernCompute(model.kern, X);
model.A = randn(model.numData, outputDim)/sqrt(model.numData+1);
model.b = randn(1, outputDim)/sqrt(model.numData+1);
