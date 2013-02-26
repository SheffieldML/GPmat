function f = lnDiffCumGaussian(u, uprime)

% LNDIFFCUMGAUSSIAN Log of the difference between two cumulative Gaussians.
%
%	Description:
%
%	F = LNDIFFCUMGAUSSIAN(U1, U2) computes the logarithm of the
%	difference between two cumulative Gaussian distributions.
%	 Returns:
%	  F - where f = log(cumGaussian(u1) - cumGaussian(u2)).
%	 Arguments:
%	  U1 - the argument of the first (positive) cumulative Gaussian.
%	  U2 - the argument of the second (negative) cumulative Gaussian.
%	
%
%	See also
%	CUMGAUSSIAN, GAUSSOVERDIFFCUMGAUSSIAN, LNCUMGAUSSIAN


%	Copyright (c) 2005, 2006 Neil D. Lawrence


f = log(gaussOverDiffCumGaussian(u, uprime, 1)+1e-300) ...
    + .5*u.*u + .5*log(2*pi);
f=-f;