function model = linearExpandParam(model, params);

% LINEAREXPANDPARAM Update linear model with vector of parameters.

startVal = 1;
endVal = model.inputDim*model.outputDim;
model.W = reshape(params(1:endVal), model.inputDim, model.outputDim);
model.b = params(endVal+1:end);