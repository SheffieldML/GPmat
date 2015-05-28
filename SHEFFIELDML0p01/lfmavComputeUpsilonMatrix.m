function [upsilonav, upsilonvp, upsilon] = lfmavComputeUpsilonMatrix(gamma, sigma2, t1, t2, mode)

% LFMAVCOMPUTEUPSILONMATRIX Upsilon matrix acce. vel. with t1, t2 limits
%
%	Description:
%
%	UPSILON = LFMAVCOMPUTEUPSILONMATRIX(GAMMA, SIGMA2, T1, T2, MODE)
%	computes a portion of the LFM kernel.
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


%	Copyright (c) 2010 Mauricio Alvarez


sigma = sqrt(sigma2);
gridt1 = repmat(t1, 1, length(t2));
gridt2 = repmat(t2', length(t1), 1);
timeGrid = gridt1 - gridt2;

if mode==0
    if nargout > 1
        [upsilonvp, upsilon] = lfmvpComputeUpsilonMatrix(gamma, sigma2, t1, t2,1);         
        upsilonav = gamma^2*upsilonvp - (4/(sqrt(pi)*sigma^3))*exp(-(timeGrid.^2)./sigma2).*...
            (timeGrid.*(gamma + 2*timeGrid/sigma2)-1);        
    else
        upsilonav = gamma^2*lfmvpComputeUpsilonMatrix(gamma, sigma2, t1, t2,1) ...
            - (4/(sqrt(pi)*sigma^3))*exp(-(timeGrid.^2)./sigma2).*...
            (timeGrid.*(gamma + 2*timeGrid/sigma2)-1);
    end    
else
    if nargout > 1
        [upsilonvp, upsilon] = lfmvpComputeUpsilonMatrix(gamma, sigma2, t1, t2, 0); 
        upsilonav = gamma^2*upsilonvp + (4/(sqrt(pi)*sigma^3))*exp(-(timeGrid.^2)./sigma2).*...
            (timeGrid.*(gamma + 2*timeGrid/sigma2)-1) ...
            - (2*gamma/(sqrt(pi)*sigma))*exp(-gamma*t1)*((gamma - 2*t2/sigma2).*exp(-(t2.^2)/sigma2)).';
    else
        upsilonav = gamma^2*lfmvpComputeUpsilonMatrix(gamma, sigma2, t1, t2, 0) ...
            + (4/(sqrt(pi)*sigma^3))*exp(-(timeGrid.^2)./sigma2).*...
            (timeGrid.*(gamma + 2*timeGrid/sigma2)-1) ...
            - (2*gamma/(sqrt(pi)*sigma))*exp(-gamma*t1)*((gamma - 2*t2/sigma2).*exp(-(t2.^2)/sigma2)).';
    end
end