function model = ppcaCreate(inputDim, outputDim, Y, options)

% PPCACREATE Density network model.
%
%	Description:
%
%	MODEL = PPCACREATE(INPUTDIMENSION, OUTPUTDIM, Y, OPTIONS) creates a
%	structure for a density network.
%	 Returns:
%	  MODEL - model structure containing the neural network specified.
%	 Arguments:
%	  INPUTDIMENSION - dimension of input data.
%	  OUTPUTDIM - dimension of target data.
%	  Y - the data to be modelled in design matrix format (as many rows
%	   as there are data points).
%	  OPTIONS - options structure as returned by ppcaCreate.
%	
%
%	See also
%	PPCAOPTIONS, MLPCREATE, RBFCREATE, KBRCREATE


%	Copyright (c) 2008 Neil D. Lawrence



model.type = 'ppca';

if size(Y, 2) ~= outputDim
  error(['Input matrix Y does not have dimension ' num2str(d)]);
end
model.X = ppcaEmbed(Y, inputDim);
model.y = Y;
model.q = inputDim;
model.d = outputDim;
model.N = size(Y, 1);
model.b = mean(model.y);
model.W = pdinv(model.X'*model.X)*model.X'*(model.y-repmat(model.b, model.N, ...
                                                  1));
yDiff = model.y - (model.X*model.W + repmat(model.b, model.N, 1));
yDiff = yDiff.*yDiff;
model.beta = model.N*model.d/sum(sum(yDiff));
