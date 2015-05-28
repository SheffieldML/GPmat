function dUpsilon = lfmvvGradientUpsilonMatrix(gamma, sigma2, t1, ...
    t2, mode, upsilon)

% LFMVVGRADIENTUPSILONMATRIX Gradient upsilon matrix vel. vel.
%
%	Description:
%
%	UPSILON = LFMVVGRADIENTUPSILONMATRIX(GAMMA, SIGMA2, T1, T2, UPSILON,
%	MODE) computes the gradient of a portion of the LFMVV kernel.
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
%	LFMVVCOMPUTEUPSILONMATRIX.M


%	Copyright (c) 2010 Mauricio Alvarez


sigma = sqrt(sigma2);

if nargin<6
    upsilon = lfmvpComputeUpsilonMatrix(gamma, sigma2, t1, t2, mode);
    if nargin <5
        mode =0;
    end
end

if mode ==0
    dUpsilon = upsilon + gamma*lfmvpGradientUpsilonMatrix(gamma, sigma2, t1, t2, mode) ...
        - (2/(sqrt(pi)*sigma))*((1-gamma*t1).*exp(-gamma*t1))*(exp(-(t2.^2)/sigma2)).';
else
    dUpsilon = -upsilon - gamma*lfmvpGradientUpsilonMatrix(gamma, sigma2, t1, t2, mode);
end

