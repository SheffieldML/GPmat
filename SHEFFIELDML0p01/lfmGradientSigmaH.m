function    g = lfmGradientSigmaH(gamma1, gamma2, sigma2, t1, t2, mode)

% LFMGRADIENTSIGMAH Gradient of the function h_i(z) with respect \sigma.
%
%	Description:
%
%	G = LFMGRADIENTSIGMAH(GAMMA1, GAMMA2, SIGMA2, T1, T2, MODE) Computes
%	the gradient of the function h_i(z) with respect to the length-scale
%	of the input "force", \sigma.
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


%	Copyright (c) 2007, 2008, Mauricio Alvarez, 2008 David Luengo


%	With modifications by Mauricio Alvarez 2008

% SEEALSO : lfmKernGradient, lfmXlfmKernGradient, lfmGradientSigmaUpsilon



% Gradient


if mode==1,
    % t1 is really t1 and t2 is really t2
    Tt1 = repmat(t1, 1, size(t2, 1));
    g = (lfmGradientSigmaUpsilon(gamma1,sigma2,t2,t1,2) - exp(-gamma2*Tt1) ...
        .*lfmGradientSigmaUpsilon(gamma1,sigma2,t2,zeros(size(t1,1)),4))/(gamma1+gamma2);
else
    % t1 is really t2 and t2 is really t1
    Tt1 = repmat(t2', size(t1, 1), 1);
    g = (lfmGradientSigmaUpsilon(gamma1,sigma2,t2,t1,1) - exp(-gamma2*Tt1) ...
        .*lfmGradientSigmaUpsilon(gamma1,sigma2,t2,zeros(size(t1,1)),3))/(gamma1+gamma2);
end

