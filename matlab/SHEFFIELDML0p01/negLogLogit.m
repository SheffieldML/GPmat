function y = negLogLogit(x)

% NEGLOGLOGIT Function which returns the negative log of the logistic function.
%
%	Description:
%
%	Y = NEGLOGLOGIT(X) computes the negative log of the logistic
%	(sigmoid) function, which is also the integral of the sigmoid
%	function.
%	 Returns:
%	  Y - the negative log of the logistic.
%	 Arguments:
%	  X - input locations.
%	
%
%	See also
%	SIGMOID


%	Copyright (c) 2006 Neil D. Lawrence


y = log(1+exp(x));