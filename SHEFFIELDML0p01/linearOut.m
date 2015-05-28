function Y = linearOut(model, X);

% LINEAROUT Obtain the output of the linear model.
%
%	Description:
%
%	Y = LINEAROUT(MODEL, X) gives the output of a linear model.
%	 Returns:
%	  Y - the output.
%	 Arguments:
%	  MODEL - the model for which the output is required.
%	  X - the input data for which the output is required.
%	
%
%	See also
%	MODELOUT, LINEARCREATE


%	Copyright (c) 2006, 2007 Neil D. Lawrence


numData = size(X, 1);
Y = X*model.W + ones(numData, 1)*model.b;
