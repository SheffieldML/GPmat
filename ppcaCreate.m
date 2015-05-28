function model = ppcaCreate(inputDim, outputDim, Y, options)

% PPCACREATE Density network model.
% FORMAT
% DESC creates a structure for a density network.
% ARG inputDimension : dimension of input data.
% ARG outputDim : dimension of target data.
% ARG Y : the data to be modelled in design matrix format (as many
% rows as there are data points).
% ARG options : options structure as returned by ppcaCreate.
% RETURN model : model structure containing the neural network
% specified.
% 
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : ppcaOptions, mlpCreate, rbfCreate, kbrCreate


% MLTOOLS

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
