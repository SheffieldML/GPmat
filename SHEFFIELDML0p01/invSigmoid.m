function y=invSigmoid(x)

% INVSIGMOID The inverse of the sigmoid function.
%
%	Description:
%
%	INVSIGMOID(X, Y) returns the inverse of the sigmoid function (which
%	takes the form y=log(x/(1-x)).
%	 Arguments:
%	  X - the input to the inverse of the sigmoid (should be between 0
%	   and 1).
%	  Y - the inverse of the sigmoid.
%	
%
%	See also
%	SIGMOID


%	Copyright (c) 2004 Neil D. Lawrence


y = log(x./(1-x));