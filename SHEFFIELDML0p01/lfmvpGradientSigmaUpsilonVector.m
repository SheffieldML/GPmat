function dUpsilonS = lfmvpGradientSigmaUpsilonVector(gamma, sigma2, t)

% LFMVPGRADIENTSIGMAUPSILONVECTOR Gradient of upsilon vector vp wrt sigma
%
%	Description:
%
%	UPSILON = LFMVPGRADIENTSIGMAUPSILONVECTOR(GAMMA, SIGMA2, T) computes
%	the gradient of a portion of the LFMVP kernel.
%	 Returns:
%	  UPSILON - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA - Gamma value for system.
%	  SIGMA2 - length scale of latent process.
%	  T - first time input (number of time points x 1).
%	
%
%	See also
%	LFMVPCOMPUTEUPSILONMATRIX.M


%	Copyright (c) 2010 Mauricio Alvarez



dUpsilon = lfmGradientSigmaUpsilonVector(gamma, sigma2, t);

dUpsilonS = -gamma*dUpsilon-(2/(sqrt(pi)*sigma2))*(1-2*(t.^2)/sigma2).* ...
    exp(-(t.^2)./sigma2);
