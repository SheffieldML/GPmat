function  g =  lfmGradientSigmaH4VV(gamma1_p, gamma1_m, sigma2, t1, ...
    preFactor, preExp)

% LFMGRADIENTSIGMAH4VV Gradient of the function h_i(z) with respect \sigma.
%
%	Description:
%
%	G = LFMGRADIENTSIGMAH4VV(GAMMA1, GAMMA2, SIGMA2, T1) Computes the
%	gradient of the function h_i(z) with respect to the length-scale of
%	the input "force", \sigma.
%	 Returns:
%	  G - Gradient of the function with respect to \sigma.
%	 Arguments:
%	  GAMMA1 - Gamma value for first system.
%	  GAMMA2 - Gamma value for second system.
%	  SIGMA2 - length scale of latent process.
%	  T1 - first time input (number of time points x 1).


%	Copyright (c) 2010 Mauricio Alvarez


g = lfmvpGradientSigmaUpsilonVector(gamma1_p,sigma2, t1)*( preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1)).' ...
    + lfmvpGradientSigmaUpsilonVector(gamma1_m,sigma2, t1)*( preExp(:,1)/preFactor(4) - preExp(:,2)/preFactor(3)).';
