function upsilon = lfmjpComputeUpsilonVector(gamma, sigma2, t)

% LFMJPCOMPUTEUPSILONVECTOR Upsilon vector jolt. pos. with t1, t2 limits
%
%	Description:
%
%	UPSILON = LFMJPCOMPUTEUPSILONVECTOR(GAMMA, SIGMA2, T1, T2) computes
%	a portion of the LFM kernel.
%	 Returns:
%	  UPSILON - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA - Gamma value for system.
%	  SIGMA2 - length scale of latent process.
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1). derivative wrt
%	   t1, mode 1 derivative wrt t2)


%	Copyright (c) 2010 Mauricio Alvarez


sigma = sqrt(sigma2);

upsilon = gamma^2*lfmvpComputeUpsilonVector(gamma, sigma2, t) ...
    + (4/(sqrt(pi)*sigma^3))*exp(-(t.^2)./sigma2).*(t.*(gamma + 2*t/sigma2) - 1);
