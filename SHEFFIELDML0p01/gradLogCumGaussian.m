function y = gradLogCumGaussian(x)

% GRADLOGCUMGAUSSIAN Gradient of the log of the cumulative Gaussian.
%
%	Description:
%
%	GRADLOGCUMGAUSSIAN(X, Y) returns the gradient of the logarithm of
%	the cumulative Gaussian. Theoretically this is simply
%	ngaussian(x)/cumGaussian(x) but there are problems in the tails of
%	the distribution. This function attempts to deal with these problems
%	without numerical error creeping in.
%	 Arguments:
%	  X - the input to the function.
%	  Y - the gradient of the log cumulative Gaussian.
%	
%
%	See also
%	NGAUSSIAN, CUMGAUSSIAN, ERFCX


%	Copyright (c) 2005 Neil D. Lawrence



y = zeros(size(x));
index = find(x>0);
y(index) = ngaussian(x(index))./cumGaussian(x(index));
index = find(x<=0);
y(index) = 1./(sqrt(2*pi)*0.5*erfcx(-sqrt(2)/2*x(index)));
