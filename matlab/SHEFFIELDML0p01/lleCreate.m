function model = lleCreate(inputDim, outputDim, Y, options)

% LLECREATE Locally linear embedding model.
%
%	Description:
%
%	MODEL = LLECREATE(LATENTDIMENSION, OUTPUTDIM, Y, OPTIONS) creates a
%	structure for a locally linear embedding.
%	 Returns:
%	  MODEL - model structure containing LLE model.
%	 Arguments:
%	  LATENTDIMENSION - dimension of latent space.
%	  OUTPUTDIM - dimension of data.
%	  Y - the data to be modelled in design matrix format (as many rows
%	   as there are data points).
%	  OPTIONS - options structure as returned by lleOptions.
%	
%
%	See also
%	LLEOPTIONS, MODELCREATE


%	Copyright (c) 2008, 2009 Neil D. Lawrence



model.type = 'lle';

if size(Y, 2) ~= outputDim
  error(['Input matrix Y does not have dimension ' num2str(d)]);
end
model.isNormalised = options.isNormalised;
model.regulariser = options.regulariser;
model.acyclic = options.acyclic;
model.k = options.numNeighbours;
model.Y = Y;
model.d = outputDim;
model.q = inputDim;
model.N = size(Y, 1);

if isfield(model, 'acyclic') && model.acyclic
  model.indices = findAcyclicNeighbours(model.Y, model.k);
else
  model.indices = findNeighbours(model.Y, model.k);
end
model.W = spalloc(model.N, model.N, model.N*model.k);
