function dUpsilon = lfmapGradientUpsilonMatrix(gamma, sigma2, t1, ...
    t2, mode, upsilon)

% LFMAPGRADIENTUPSILONMATRIX Gradient upsilon matrix accel. pos.
%
%	Description:
%
%	UPSILON = LFMAPGRADIENTUPSILONMATRIX(GAMMA, SIGMA2, T1, T2, UPSILON,
%	MODE) computes the gradient of a portion of the LFMAP kernel.
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
%	LFMAPCOMPUTEUPSILONMATRIX.M


%	Copyright (c) 2010 Mauricio Alvarez


sigma = sqrt(sigma2);
gridt1 = repmat(t1, 1, length(t2));
gridt2 = repmat(t2', length(t1), 1);
timeGrid = gridt1 - gridt2;

if nargin<6
    upsilon = lfmComputeUpsilonMatrix(gamma, sigma2, t1, t2);
    if nargin <5
        mode =0;
    end
end

if mode ==0
    dUpsilon = 2*gamma*upsilon + gamma^2*lfmGradientUpsilonMatrix(gamma, sigma2, t1, t2) ...
        - (2/(sqrt(pi)*sigma))*exp(-(timeGrid.^2)./sigma2);
else
    dUpsilon = 2*gamma*upsilon + gamma^2*lfmGradientUpsilonMatrix(gamma, sigma2, t1, t2) ...
        - (2/(sqrt(pi)*sigma))*exp(-(timeGrid.^2)./sigma2) ...
        - (2/(sqrt(pi)*sigma))*(t1.*exp(-gamma*t1))*((gamma - 2*t2/sigma2).*exp(-t2.^2/sigma2)).' ...
        + (2/(sqrt(pi)*sigma))*exp(-gamma*t1)*(exp(-t2.^2/sigma2)).';
end

