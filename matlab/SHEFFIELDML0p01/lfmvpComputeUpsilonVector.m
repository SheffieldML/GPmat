function [upsilonvp, upsilon] = lfmvpComputeUpsilonVector(gamma, sigma2, t1 ,mode)

% LFMVPCOMPUTEUPSILONVECTOR Upsilon vector for vel. pos. with t1 limit
%
%	Description:
%
%	UPSILON = LFMVPCOMPUTEUPSILONVECTOR(GAMMA, SIGMA2, T1, MODE)
%	computes a portion of the LFMVP kernel.
%	 Returns:
%	  UPSILON - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA - Gamma value for system.
%	  SIGMA2 - length scale of latent process.
%	  T1 - first time input (number of time points x 1).
%	  MODE - operation mode, according to the derivative (mode 0,
%	   derivative wrt t1, mode 1 derivative wrt t2)
%	
%
%	See also
%	LFMVPCOMPUTEUPSILONMATRIX.M


%	Copyright (c) 2010 Mauricio A. Alvarez


if nargin<4
    mode = 0;
end

sigma = sqrt(sigma2);
 
if mode==0 
    if nargout > 1
        upsilon = lfmComputeUpsilonVector(gamma, sigma2, t1); 
        upsilonvp = -gamma*upsilon + (2/(sqrt(pi)*sigma))*exp(-(t1.^2)/sigma2);
    else
        upsilonvp = -gamma*lfmComputeUpsilonVector(gamma, sigma2, t1) ...
            + (2/(sqrt(pi)*sigma))*exp(-(t1.^2)/sigma2);
    end
else
end