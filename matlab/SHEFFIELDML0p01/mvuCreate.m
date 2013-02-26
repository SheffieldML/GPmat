function model = mvuCreate(inputDim, outputDim, Y, options)

% MVUCREATE Maximum variance unfolding embedding model.
%
%	Description:
%
%	MODEL = MVUCREATE(LATENTDIMENSION, OUTPUTDIM, Y, OPTIONS) creates a
%	structure for an mvu model.
%	 Returns:
%	  MODEL - model structure containing the MVU model.
%	 Arguments:
%	  LATENTDIMENSION - dimension of latent space.
%	  OUTPUTDIM - dimension of data.
%	  Y - the data to be modelled in design matrix format (as many rows
%	   as there are data points).
%	  OPTIONS - options structure as returned by mvuOptions.
%	
%
%	See also
%	MVUOPTIONS, LLECREATE, ISOMAPCREATE, MODELCREATE


%	Copyright (c) 2009 Neil D. Lawrence



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