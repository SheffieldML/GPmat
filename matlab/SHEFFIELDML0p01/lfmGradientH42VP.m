function    g = lfmGradientH42VP(preFactor, preFactorGrad, gradThetaGamma, ...
    preExp, preExpg, preExpt, compUpsilon1, compUpsilon2)

% LFMGRADIENTH42VP Gradient of the function h_i(z) with respect to some of the
%
%	Description:
%	hyperparameters of the kernel: m_k, C_k, D_k, m_r, C_r or D_r.
%
%	G = LFMGRADIENTH42VP(GAMMA1, GAMMA2, SIGMA2, GRADTHETAGAMMA, T1, T2)
%	Computes the gradient of the function h_i(z) with respect to some of
%	the parameters of the system (mass, spring or damper).
%	 Returns:
%	  G - Gradient of the function with respect to the desired
%	   parameter.
%	 Arguments:
%	  GAMMA1 - Gamma value for first system.
%	  GAMMA2 - Gamma value for second system.
%	  SIGMA2 - length scale of latent process.
%	  GRADTHETAGAMMA - Vector with the gradient of gamma1 and gamma2
%	   with respect to the desired parameter.
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1)
%	
%
%	See also
%	LFMKERNGRADIENT, LFMXLFMKERNGRADIENT, LFMGRADIENTUPSILON


%	Copyright (c) 2010 Mauricio A. Alvarez



% Gradient evaluation


g = compUpsilon1*(( - preExpg(:,2)/preFactorGrad(3) - preExpt(:,2)/preFactor(3) + preExp(:,2)/preFactor(3))*gradThetaGamma(2)...
    + ( preExpg(:,1)/preFactorGrad(1) + preExpt(:,1)/preFactor(1) - preExp(:,1)/preFactor(1))*gradThetaGamma(1)).'...
    - compUpsilon2*(( - preExpg(:,2)/preFactorGrad(4) - preExpt(:,2)/preFactor(4) + preExp(:,2)/preFactor(4))*gradThetaGamma(2)...
    + ( preExpg(:,1)/preFactorGrad(2) + preExpt(:,1)/preFactor(2) - preExp(:,1)/preFactor(2))*gradThetaGamma(1)).';

    

