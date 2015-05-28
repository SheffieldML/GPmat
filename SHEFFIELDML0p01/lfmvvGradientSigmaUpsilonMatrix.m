function dUpsilonS = lfmvvGradientSigmaUpsilonMatrix(gamma, sigma2, ...
    t1, t2, mode)

% LFMVVGRADIENTSIGMAUPSILONMATRIX Gradient of upsilon matrix vv wrt sigma
%
%	Description:
%
%	UPSILON = LFMVVGRADIENTSIGMAUPSILONMATRIX(GAMMA, SIGMA2, T1, T2,
%	MODE) computes the gradient wrt sigma of a portion of the LFMVV
%	kernel.
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
%	LFMVVCOMPUTEUPSILONMATRIX.M


%	Copyright (c) 2010 Mauricio Alvarez


gridt1 = repmat(t1, 1, length(t2));
gridt2 = repmat(t2', length(t1), 1);
timeGrid = gridt1 - gridt2;

dUpsilon = lfmvpGradientSigmaUpsilonMatrix(gamma, sigma2, t1, t2, mode);

if mode == 0
    dUpsilonS = gamma*dUpsilon - 4*timeGrid/(sqrt(pi)*sigma2^2).* ...
        exp(-(timeGrid.^2)./sigma2).*(3 - 2*(timeGrid.^2)/sigma2) ...
        + 2*gamma/(sqrt(pi)*sigma2)*exp(-gamma*t1)*((1-2*t2.^2/sigma2).* ...
        exp(-(t2.^2)/sigma2)).';
else
    dUpsilonS = -gamma*dUpsilon - 4*timeGrid/(sqrt(pi)*sigma2^2).* ...
        exp(-(timeGrid.^2)./sigma2).*(3 - 2*(timeGrid.^2)/sigma2);
end
