function y = gaussSamp(Sigma, numSamps)

% GAUSSSAMP Sample from a Gaussian with a given covariance.
%
%	Description:
%
%	Y = GAUSSSAMP(SIGMA, NUMSAMPS) samples a given number of samples
%	from a Gaussian with a given covariance matrix.
%	 Returns:
%	  Y - the samples from the Gaussian
%	 Arguments:
%	  SIGMA - the covariance of the Gaussian to sample from.
%	  NUMSAMPS - the number of samples to take from Gaussian.
%	
%
%	See also
%	RANDN, EIG


%	Copyright (c) 2005 Neil D. Lawrence


[U, V] = eig(Sigma);
dims = size(Sigma, 1);
y = randn(numSamps, dims);
y = y*diag(sqrt(diag(V)));
y = y*U';