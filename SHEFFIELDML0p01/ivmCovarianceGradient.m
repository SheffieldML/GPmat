function g = ivmCovarianceGradient(invK, m)

% IVMCOVARIANCEGRADIENT The gradient of the likelihood approximation wrt the covariance.
%
%	Description:
%
%	G = IVMCOVARIANCEGRADIENT(INVK, M) returns the gradient of the log
%	likelihood approximation with respect to the covariance function
%	(i.e. kernel) evaluated at the active points.
%	 Returns:
%	  G - the gradient of the log likelihood approximation with respect
%	   to the covariance function elements.
%	 Arguments:
%	  INVK - the inverse of the covariance function evaluated at the
%	   active points.
%	  M - the difference between the target and the mean value.
%	
%
%	See also
%	IVMCREATE


%	Copyright (c) 2004, 2005 Neil D. Lawrence


invKm = invK*m;

g = -invK + invKm*invKm';
g= g*.5;
