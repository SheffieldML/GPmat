function model = linearParamInit(model)

% LINEARPARAMINIT Initialise the parameters of an LINEAR model.
%
%	Description:
%
%	MODEL = LINEARPARAMINIT(MODEL) sets the initial weight vectors and
%	biases to small random values.
%	 Returns:
%	  MODEL - the initialised model.
%	 Arguments:
%	  MODEL - the input model to initialise.
%	
%
%	See also
%	MODELPARAMINIT, LINEARCREATE


%	Copyright (c) 2006 Neil D. Lawrence


model.W = randn(model.inputDim, model.outputDim)/sqrt(model.inputDim + 1);
model.b = randn(1, model.outputDim)/sqrt(model.inputDim + 1);
model.beta = 1;