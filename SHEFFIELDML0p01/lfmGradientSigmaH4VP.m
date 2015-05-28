function    g = lfmGradientSigmaH4VP(gamma1, gamma2, sigma2, t1, ...
    preFactor, preExp, mode)

% LFMGRADIENTSIGMAH4VP Gradient of the function h_i(z) with respect \sigma.
%
%	Description:
%
%	G = LFMGRADIENTSIGMAH4VP(GAMMA1, GAMMA2, SIGMA2, T1, T2, MODE)
%	Computes the gradient of the function h_i(z) with respect to the
%	length-scale of the input "force", \sigma.
%	 Returns:
%	  G - Gradient of the function with respect to \sigma.
%	 Arguments:
%	  GAMMA1 - Gamma value for first system.
%	  GAMMA2 - Gamma value for second system.
%	  SIGMA2 - length scale of latent process.
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1).
%	  MODE - indicates in which way the vectors t1 and t2 must be
%	   transposed
%	
%
%	See also
%	LFMKERNGRADIENT, LFMXLFMKERNGRADIENT, LFMGRADIENTSIGMAUPSILON


%	Copyright (c) 2010 Mauricio Alvarez


if mode==0
    g =  lfmvpGradientSigmaUpsilonVector(gamma1,sigma2, t1)*( preExp(:,1)/preFactor(1) - preExp(:,2)/preFactor(2)).' ...
        + lfmvpGradientSigmaUpsilonVector(gamma2,sigma2, t1)*( preExp(:,2)/preFactor(3) - preExp(:,1)/preFactor(4)).';
else
    g =  lfmGradientSigmaUpsilonVector(gamma1,sigma2, t1)*(  preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1)).' ...
        + lfmGradientSigmaUpsilonVector(gamma2,sigma2, t1)*( preExp(:,1)/preFactor(4) - preExp(:,2)/preFactor(3)).';
end



