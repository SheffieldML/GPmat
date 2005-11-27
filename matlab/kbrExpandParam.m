function model = kbrExpandParam(model, params);

% KBREXPANDPARAM Update kernel based regression model with vector of parameters.

startVal = 1;
endVal = model.numData*model.outputDim;
model.A = reshape(params(1:endVal), model.numData, model.outputDim);
model.b = params(endVal+1:end);