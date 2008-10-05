function model = lleCreate(inputDim, outputDim, Y, options)

% LLECREATE Locally linear embedding model.
% FORMAT
% DESC creates a structure for a locally linear embedding.
% ARG inputDimension : dimension of latent space.
% ARG outputDim : dimension of data.
% ARG Y : the data to be modelled in design matrix format (as many
% rows as there are data points).
% ARG options : options structure as returned by lleCreate.
% RETURN model : model structure containing the neural network
% specified.
% 
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : lleOptions, modelCreate


% MLTOOLS

model.type = 'lle';

if size(Y, 2) ~= outputDim
  error(['Input matrix Y does not have dimension ' num2str(d)]);
end
model.k = options.numNeighbours;
model.Y = Y;
model.d = outputDim;
model.q = options.latentDim;
model.N = size(Y, 1);