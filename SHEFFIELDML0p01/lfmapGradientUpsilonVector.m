function dUpsilon = lfmapGradientUpsilonVector(gamma, sigma2, t, upsilon)

% LFMAPGRADIENTUPSILONVECTOR Gradient upsilon vector accel. pos.
%
%	Description:
%
%	UPSILON = LFMAPGRADIENTUPSILONVECTOR(GAMMA, SIGMA2, T1, T2, UPSILON)
%	computes the gradient of a portion of the LFM kernel.
%	 Returns:
%	  UPSILON - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA - Gamma value for system.
%	  SIGMA2 - length scale of latent process.
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1).
%	  UPSILON - precomputation of the upsilon matrix.
%	
%
%	See also
%	LFMAPCOMPUTEUPSILONVECTOR.M


%	Copyright (c) 2010 Mauricio Alvarez


sigma = sqrt(sigma2);

if nargin<4
    upsilon = lfmComputeUpsilonVector(gamma, sigma2, t);
end

dUpsilon = 2*gamma*upsilon + gamma^2*lfmGradientUpsilonVector(gamma, sigma2, t) ...
        - (2/(sqrt(pi)*sigma))*exp(-(t.^2)./sigma2);

