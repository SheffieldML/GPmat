function    g = lfmGradientSigmaH4(gamma1, gamma2, sigma2, t1, ...
    preFactor, preExp, mode, term)

% LFMGRADIENTSIGMAH4 Gradient of the function h_i(z) with respect \sigma.
%
%	Description:
%
%	G = LFMGRADIENTSIGMAH4(GAMMA1, GAMMA2, SIGMA2, T1, T2, MODE)
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
%
%	See also
%	LFMKERNGRADIENT, LFMXLFMKERNGRADIENT, LFMGRADIENTSIGMAUPSILON


%	Copyright (c) 2007, 2008 David Luengo
%	Copyright (c) 2008 Mauricio Alvarez



% Gradient

if nargin<8
    term =[];
end

if ~mode
    if ~term
        g = lfmGradientSigmaUpsilonVector(gamma1,sigma2, t1)*( preExp/preFactor(1) - conj(preExp)/preFactor(2)).';
    else
        gradupsilon = lfmGradientSigmaUpsilonVector(gamma1,sigma2, t1);
        g = gradupsilon*(preExp/preFactor(1)).' - conj(gradupsilon)*(preExp/preFactor(2)).';
    end
else
    g =  lfmGradientSigmaUpsilonVector(gamma1,sigma2,t1)*( preExp(:,1)/preFactor(1) - preExp(:,2)/preFactor(2)).' ...
        + lfmGradientSigmaUpsilonVector(gamma2,sigma2,t1)*( preExp(:,2)/preFactor(3) - preExp(:,1)/preFactor(4)).';
end





