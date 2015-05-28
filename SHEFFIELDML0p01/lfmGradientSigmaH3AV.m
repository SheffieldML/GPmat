function g =  lfmGradientSigmaH3AV(gamma1_p, gamma1_m, sigma2, t1, ...
    t2, preFactor, mode)

% LFMGRADIENTSIGMAH3AV Gradient of the function h_i(z) with respect \sigma.
%
%	Description:
%
%	G = LFMGRADIENTSIGMAH3AV(GAMMA1, GAMMA2, SIGMA2, T1, T2, MODE)
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


%	Copyright (c) 2010 Mauricio Alvarez


g = preFactor(1)*lfmavGradientSigmaUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode) ...
    + preFactor(2)*lfmavGradientSigmaUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);