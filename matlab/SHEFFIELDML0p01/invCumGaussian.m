function y = invCumGaussian(x)

% INVCUMGAUSSIAN Computes inverse of the cumulative Gaussian.
%
%	Description:
%
%	Y = INVCUMGAUSSIAN(X) computes the inverse of the cumulative
%	Gaussian.
%	 Returns:
%	  Y - the inverse of the cumulative Gaussian.
%	 Arguments:
%	  X - value between 0 and 1 to map onto the real line.
%	
%
%	See also
%	CUMGAUSSIAN, ERFINV


%	Copyright (c) 2005 Neil D. Lawrence


y = erfinv(x*2 - 1)*2/sqrt(2);