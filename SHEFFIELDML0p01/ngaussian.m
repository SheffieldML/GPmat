function y = ngaussian(x)

% NGAUSSIAN Compute a Gaussian with mean 0 and variance 1.
%
%	Description:
%
%	Y = NGAUSSIAN(X) computes a the likelihood of a normalised Gaussian
%	distribution, i.e. with mean 0 and variance 1.
%	 Returns:
%	  Y - probability of the input values under the Gaussian.
%	 Arguments:
%	  X - input value(s) for which to compute the distribution.
%	
%
%	See also
%	CUMGAUSSIAN, GAUSSOVERDIFFCUMGAUSSIAN


%	Copyright (c) 2004 Neil D. Lawrence


x2 = x.*x;
y = exp(-.5*x2);
y = y/sqrt(2*pi);
