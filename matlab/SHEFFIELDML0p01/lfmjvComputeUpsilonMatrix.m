function upsilonjv = lfmjvComputeUpsilonMatrix(gamma, sigma2, t1, t2, mode)

% LFMJVCOMPUTEUPSILONMATRIX Upsilon matrix jolt. vel. with t1, t2 limits
%
%	Description:
%
%	UPSILON = LFMJVCOMPUTEUPSILONMATRIX(GAMMA, SIGMA2, T1, T2, MODE)
%	computes a portion of the LFMJV kernel.
%	 Returns:
%	  UPSILON - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA - Gamma value for system.
%	  SIGMA2 - length scale of latent process.
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1).
%	  MODE - operation mode, according to the derivative (mode 0,
%	   derivative wrt t1, mode 1 derivative wrt t2)
%	
%
%	See also
%	LFMCOMPUTEUPSILONMATRIX.F, LFMVPCOMPUTEUPSILONMATRIX.M


%	Copyright (c) 2010 Mauricio Alvarez


sigma = sqrt(sigma2);
gridt1 = repmat(t1, 1, length(t2));
gridt2 = repmat(t2', length(t1), 1);
timeGrid = gridt1 - gridt2;

if mode==0
    upsilonjv = gamma^2*lfmvvComputeUpsilonMatrix(gamma, sigma2, t1, t2, mode) ...
        - (4/(sqrt(pi)*sigma^3))*exp(-(timeGrid.^2)./sigma2).* ...
        ((gamma + (2*timeGrid)/sigma2).*(1-(2*timeGrid.^2)/sigma2) + 4*timeGrid/sigma2);
else
    upsilonjv = gamma^2*lfmvvComputeUpsilonMatrix(gamma, sigma2, t1, t2, mode) ...
        - (4/(sqrt(pi)*sigma^3))*exp(-(timeGrid.^2)./sigma2).* ...
        ((gamma + (2*timeGrid)/sigma2).*(1-(2*timeGrid.^2)/sigma2) + 4*timeGrid/sigma2) ...
        + ((4*gamma)/(sqrt(pi)*sigma^3))*exp(-gamma*t1)*...
        ((t2.*(gamma - 2*t2/sigma2) + 1).*exp(-t2.^2/sigma2)).';
end