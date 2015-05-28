function model = kbrOptimise(model, X, Y, varargin)

% KBROPTIMISE Optimise a KBR model.
%
%	Description:
%
%	MODEL = KBROPTIMISE(MODEL, X, Y) optimises a kernel based regression
%	model using a least squares fit.
%	 Returns:
%	  MODEL - the optimised model.
%	 Arguments:
%	  MODEL - the model to be optimised.
%	  X - the input data locations for the optimisation.
%	  Y - the target data locations for the optimisation.
%	
%
%	See also
%	KBRCREATE, MODELOPTIMISE


%	Copyright (c) 2005, 2006 Neil D. Lawrence



model.numData = size(X, 1);
model.K = kernCompute(model.kern, X);
model.X = X;

model.bias = mean(Y, 1);
model.A = pdinv(model.K)*(Y-repmat(model.bias, model.numData, 1));
