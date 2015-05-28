function y = cumGamma(x, a, b)

% CUMGAMMA Cumulative distribution for gamma.
%
%	Description:
%
%	P = CUMGAMMA(X) computes the cumulative gamma distribution.
%	 Returns:
%	  P - output probability.
%	 Arguments:
%	  X - input value.
%	
%
%	See also
%	GAMMAINC, GAMMA


%	Copyright (c) 2008 Neil D. Lawrence


y = gammainc(x*b, a);
