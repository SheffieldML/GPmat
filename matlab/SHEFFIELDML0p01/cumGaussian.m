function y = cumGaussian(x)

% CUMGAUSSIAN Cumulative distribution for Gaussian.
%
%	Description:
%
%	P = CUMGAUSSIAN(X) computes the cumulative Gaussian distribution.
%	 Returns:
%	  P - output probability.
%	 Arguments:
%	  X - input value.
%	
%
%	See also
%	LNCUMGAUSSIAN, LNDIFFCUMGAUSSIAN, ERF


%	Copyright (c) 2004 Neil D. Lawrence


y = 0.5*(1+erf(sqrt(2)/2*x));
