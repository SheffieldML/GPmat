function upsilon = lfmvvComputeUpsilonDiagVector(gamma, sigma2, t, mode)

% LFMVVCOMPUTEUPSILONDIAGVECTOR Upsilon vector vel. vel. with t1 = t2
%
%	Description:
%
%	UPSILON = LFMVVCOMPUTEUPSILONDIAGVECTOR(GAMMA, SIGMA2, T, MODE)
%	computes a portion of the LFMVV kernel.
%	 Returns:
%	  UPSILON - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA - Gamma value for system.
%	  SIGMA2 - length scale of latent process.
%	  T - first time input (number of time points x 1).
%	  MODE - operation mode, according to the derivative (mode 0,
%	   derivative wrt t1, mode 1 derivative wrt t2)
%	
%
%	See also
%	LFMCOMPUTEUPSILONMATRIX.F, LFMVPCOMPUTEUPSILONMATRIX.M


%	Copyright (c) 2010 Mauricio Alvarez


sigma = sqrt(sigma2);

if mode==0
    upsilon = gamma*lfmvpComputeUpsilonDiagVector(gamma, sigma2, t, mode) ...
        - (2*gamma/(sqrt(pi)*sigma))*exp(-gamma*t).*(exp(-(t.^2)/sigma2));
else
    upsilon = -gamma*lfmvpComputeUpsilonDiagVector(gamma, sigma2, t, mode);
end