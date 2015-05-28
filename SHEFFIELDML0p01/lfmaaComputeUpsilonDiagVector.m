function upsilon = lfmaaComputeUpsilonDiagVector(gamma, sigma2, t, mode)

% LFMAACOMPUTEUPSILONDIAGVECTOR Diag. of Upsilon matrix acce. accel.
%
%	Description:
%
%	UPSILON = LFMAACOMPUTEUPSILONDIAGVECTOR(GAMMA, SIGMA2, T, MODE)
%	computes a portion of the LFMAXLFMA kernel.
%	 Returns:
%	  UPSILON - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA - Gamma value for system.
%	  SIGMA2 - length scale of latent process.
%	  T - first time input (number of time points x 1).
%	  MODE - operation mode, according to the derivative (mode 0,
%	   derivative wrt t1, mode 1 derivative wrt t2)


%	Copyright (c) 2010 Mauricio Alvarez


sigma = sqrt(sigma2);

if mode==0   
    upsilon = gamma^2*lfmapComputeUpsilonDiagVector(gamma, sigma2, t, 1) ...
        + (2/(sqrt(pi)*sigma)*(gamma*2/sigma2));   
else
    upsilon = gamma^2*lfmapComputeUpsilonDiagVector(gamma, sigma2, t, 0) ...
        + (2/(sqrt(pi)*sigma)*(gamma*2/sigma2)) ...
        + (2*gamma^2/(sqrt(pi)*sigma))*exp(-gamma*t).*((gamma - 2*t/sigma2).*exp(-(t.^2)/sigma2));
end