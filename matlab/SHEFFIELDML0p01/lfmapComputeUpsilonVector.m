function [upsilonap, upsilon] = lfmapComputeUpsilonVector(gamma, sigma2, t1)

% LFMAPCOMPUTEUPSILONVECTOR Upsilon vector for acce. pos. with t1 limit
%
%	Description:
%
%	UPSILON = LFMAPCOMPUTEUPSILONVECTOR(GAMMA, SIGMA2, T1) computes a
%	portion of the LFMAP kernel.
%	 Returns:
%	  UPSILON - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA - Gamma value for system.
%	  SIGMA2 - length scale of latent process.
%	  T1 - first time input (number of time points x 1).
%	
%
%	See also
%	LFMAPCOMPUTEUPSILONMATRIX.M


%	Copyright (c) 2010 Mauricio Alvarez


sigma = sqrt(sigma2);

if nargout > 1
    upsilon = lfmComputeUpsilonVector(gamma, sigma2, t1);
    upsilonap = gamma^2*upsilon - (2/(sqrt(pi)*sigma))*exp(-(t1.^2)/sigma2).*(gamma + 2*t1/sigma2);
else
    upsilonap = gamma^2*lfmComputeUpsilonVector(gamma, sigma2, t1) ...
        - (2/(sqrt(pi)*sigma))*exp(-(t1.^2)/sigma2).*(gamma + 2*t1/sigma2);
end