function model = leCreate(inputDim, outputDim, Y, options)

% LECREATE Laplacian eigenmap model.
% FORMAT
% DESC creates a structure for a Laplacian eigenmap model.
% ARG latentDimension : dimension of latent space.
% ARG outputDim : dimension of data.
% ARG Y : the data to be modelled in design matrix format (as many
% rows as there are data points).
% ARG options : options structure as returned by leOptions.
% RETURN model : model structure containing LE model.
% 
% COPYRIGHT : Neil D. Lawrence, 2009
%
% SEEALSO : leOptions, modelCreate


% MLTOOLS

model.type = 'le';

if size(Y, 2) ~= outputDim
  error(['Input matrix Y does not have dimension ' num2str(d)]);
end
model.isNormalised = options.isNormalised;
model.weightType = options.weightType;
model.weightScale = options.weightScale;
model.regulariser = options.regulariser;
model.k = options.numNeighbours;
model.Y = Y;
model.d = outputDim;
model.q = inputDim;
model.N = size(Y, 1);

[model.indices, D2] = findNeighbours(model.Y, model.k);
model.L = spalloc(model.N, model.N, model.N*model.k);
