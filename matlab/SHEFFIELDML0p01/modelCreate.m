function model = modelCreate(type, numIn, numOut, varargin)

% MODELCREATE Create a model of the specified type.
%
%	Description:
%
%	MODEL = MODELCREATE(TYPE, NUMIN, NUMOUT, P3,...) creates a model of
%	the given type.
%	 Returns:
%	  MODEL - the model created.
%	 Arguments:
%	  TYPE - the type of the model to create, for example, 'kbr' for
%	   kernel based regression, 'mlp' for multi-layer perceptron,
%	   'linear' for a linear model.
%	  NUMIN - number of inputs to the model (or latent dimensions for
%	   latent variable models.
%	  NUMOUT - number of outputs from the model (or data dimensions for
%	   latent variable models.
%	  P3,... - optional arguments to be passed to the model creation
%	   code.
%	
%
%	See also
%	MODELEXPANDPARAM, MODELEXTRACTPARAM


%	Copyright (c) 2005, 2006 Neil D. Lawrence


fhandle = str2func([type, 'Create']);
model = fhandle(numIn, numOut, varargin{:});