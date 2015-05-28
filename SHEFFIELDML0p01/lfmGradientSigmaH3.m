function    g = lfmGradientSigmaH3(gamma1, gamma2, sigma2, t1, t2, preFactor, mode, term)

% LFMGRADIENTSIGMAH3 Gradient of the function h_i(z) with respect \sigma.
%
%	Description:
%
%	G = LFMGRADIENTSIGMAH3(GAMMA1, GAMMA2, SIGMA2, T1, T2, MODE)
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


%	Copyright (c) 2007, 2008, Mauricio Alvarez, 2008 David Luengo


%	With modifications by Mauricio Alvarez 2008

% SEEALSO : lfmKernGradient, lfmXlfmKernGradient, lfmGradientSigmaUpsilon



% Gradient

if nargin<8
    term = [];
end

if ~mode
    if ~term
        g = preFactor*lfmGradientSigmaUpsilonMatrix(gamma1,sigma2, t1,t2);
    else
        gradupsilon = lfmGradientSigmaUpsilonMatrix(gamma1,sigma2, t1,t2);
        g = -preFactor(1)*gradupsilon + preFactor(2)*conj(gradupsilon);
    end
else
    g =  preFactor(1)*lfmGradientSigmaUpsilonMatrix(gamma1,sigma2,t1,t2) + ...
        preFactor(2)*lfmGradientSigmaUpsilonMatrix(gamma2,sigma2,t1,t2);
end



