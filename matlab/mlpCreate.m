function model = mlpCreate(inputDim, outputDim, hiddenDim, activeFunc)

% MLPCREATE Wrapper for NETLAB's mlp `net'.

if nargin < 5
  hiddenDim = 20;
  if nargin < 4
    activeFunc = 'linear';
  end
end

model = mlp(inputDim, hiddenDim, outputDim, activeFunc);
model.numParams = model.nwts;
model.inputDim = inputDim;
model.outputDim = outputDim;
