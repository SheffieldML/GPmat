function model = rbfCreate(inputDim, outputDim, hiddenDim, activeFunc, outFunc)

% RBFCREATE Wrapper for NETLAB's rbf `net'.

if nargin < 5
  activeFunc = 'linear';
  if nargin < 4
    activeFunc = 'gaussian';
    if nargin < 3
      hiddenDim = 20;
    end
  end
end

model = rbf(inputDim, hiddenDim, outputDim, activeFunc, outFunc);
model.numParams = model.nwts;
model.inputDim = inputDim;
model.outputDim = outputDim;
