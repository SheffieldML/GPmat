function model = mvuCreate(inputDim, outputDim, Y, options)

% MVUCREATE Maximum variance unfolding embedding model.
% FORMAT
% DESC creates a structure for an mvu model.
% ARG latentDimension : dimension of latent space.
% ARG outputDim : dimension of data.
% ARG Y : the data to be modelled in design matrix format (as many
% rows as there are data points).
% ARG options : options structure as returned by mvuOptions.
% RETURN model : model structure containing the MVU model.
% 
% COPYRIGHT : Neil D. Lawrence, 2009
%
% SEEALSO : mvuOptions, lleCreate, isomapCreate, modelCreate


% MLTOOLS

model.type = 'mvu';

if size(Y, 2) ~= outputDim
  error(['Input matrix Y does not have dimension ' num2str(d)]);
end
model.k = options.numNeighbours;
model.Y = Y;
model.solver = options.solver;
model.d = outputDim;
model.q = inputDim;
model.N = size(Y, 1);
