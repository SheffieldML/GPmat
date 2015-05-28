function dUpsilon = lfmvpGradientUpsilonVector(gamma, sigma2, t, upsilon)

% LFMVPGRADIENTUPSILONVECTOR Gradient upsilon vector vel. pos.
%
%	Description:
%
%	UPSILON = LFMVPGRADIENTUPSILONVECTOR(GAMMA, SIGMA2, T1, T2, UPSILON,
%	MODE) computes the gradient of a portion of the LFM kernel.
%	 Returns:
%	  UPSILON - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA - Gamma value for system.
%	  SIGMA2 - length scale of latent process.
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1).
%	  UPSILON - precomputation of the upsilon matrix.
%	  MODE - operation mode, according to the derivative (mode 0,
%	   derivative wrt t1, mode 1 derivative wrt t2)
%	
%
%	See also
%	LFMVPCOMPUTEUPSILONMATRIX.M


%	Copyright (c) 2010 Mauricio A. Alvarez


if nargin<4
    upsilon = lfmComputeUpsilonVector(gamma, sigma2, t);
end

dUpsilon = -upsilon - gamma*lfmGradientUpsilonVector(gamma, sigma2, t);

