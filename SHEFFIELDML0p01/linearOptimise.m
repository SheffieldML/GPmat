function model = linearOptimise(model, X, Y, varargin)

% LINEAROPTIMISE Optimise a linear model.
%
%	Description:
%
%	MODEL = LINEAROPTIMISE(MODEL, X, Y) optimises a linear model using
%	least squares.
%	 Returns:
%	  MODEL - the optimised model.
%	 Arguments:
%	  MODEL - the model to be optimised.
%	  X - the input data locations for the optimisation.
%	  Y - the target data locations for the optimisation.
%	
%
%	See also
%	LINEARCREATE, MODELOPTIMISE


%	Copyright (c) 2005, 2006, 2008 Neil D. Lawrence


N = size(X, 1);
Xo = [X ones(N, 1)];
W = pdinv(Xo'*Xo)*Xo'*Y;
model.b = W(end, :);
model.W = W(1:end-1, :);
%model.b = mean(Y); %W(end, :);
%model.W = pdinv(X'*X)*X'*(Y - repmat(model.b, size(Y, 1), 1));
centred = Y - linearOut(model, X);
model.beta = N./sum(centred.*centred, 1);
