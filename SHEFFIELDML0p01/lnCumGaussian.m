function y = lnCumGaussian(x)

% LNCUMGAUSSIAN log cumulative distribution for the normalised Gaussian.
%
%	Description:
%
%	Y = LNCUMGAUSSIAN(X) computes the logarithm of the cumulative
%	Gaussian distribution.
%	 Returns:
%	  Y - log probability of the value under the cumulative Gaussian.
%	 Arguments:
%	  X - input position.
%	
%
%	See also
%	ERF, ERFCX, CUMGAUSSIAN, LNDIFFCUMGAUSSIAN, GAUSSOVERDIFFCUMGAUSSIAN


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


index = find(x< 0);
if length(index)
  y(index) = -.5*x(index).*x(index) + log(.5) + log(erfcx(-sqrt(2)/2* ...
                                                    x(index)));
end
index = find(x>=0);
if length(index)
  y(index) = log(cumGaussian(x(index)));
end
y=reshape(y, size(x));