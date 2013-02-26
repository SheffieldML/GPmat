function model = isomapCreate(inputDim, outputDim, Y, options)

% ISOMAPCREATE isomap embedding model.
%
%	Description:
%
%	MODEL = ISOMAPCREATE(LATENTDIMENSION, OUTPUTDIM, Y, OPTIONS) creates
%	a structure for an isomap model.
%	 Returns:
%	  MODEL - model structure containing the isomap model.
%	 Arguments:
%	  LATENTDIMENSION - dimension of latent space.
%	  OUTPUTDIM - dimension of data.
%	  Y - the data to be modelled in design matrix format (as many rows
%	   as there are data points).
%	  OPTIONS - options structure as returned by isomapOptions.
%	
%
%	See also
%	ISOMAPOPTIONS, LLECREATE, MVUCREATE, MODELCREATE


%	Copyright (c) 2009 Neil D. Lawrence



model.type = 'isomap';

if size(Y, 2) ~= outputDim
  error(['Input matrix Y does not have dimension ' num2str(d)]);
end
model.k = options.numNeighbours;
model.Y = Y;
model.d = outputDim;
model.q = inputDim;
model.N = size(Y, 1);