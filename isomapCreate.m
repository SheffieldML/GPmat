function model = isomapCreate(inputDim, outputDim, Y, options)

% ISOMAPCREATE isomap embedding model.
% FORMAT
% DESC creates a structure for an isomap model.
% ARG latentDimension : dimension of latent space.
% ARG outputDim : dimension of data.
% ARG Y : the data to be modelled in design matrix format (as many
% rows as there are data points).
% ARG options : options structure as returned by isomapOptions.
% RETURN model : model structure containing the isomap model.
% 
% COPYRIGHT : Neil D. Lawrence, 2009
%
% SEEALSO : isomapOptions, lleCreate, mvuCreate, modelCreate


% MLTOOLS

model.type = 'isomap';

if size(Y, 2) ~= outputDim
  error(['Input matrix Y does not have dimension ' num2str(d)]);
end
model.k = options.numNeighbours;
model.Y = Y;
model.d = outputDim;
model.q = inputDim;
model.N = size(Y, 1);
