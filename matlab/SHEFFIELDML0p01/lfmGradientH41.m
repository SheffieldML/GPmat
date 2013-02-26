function    g = lfmGradientH41(preFactor, preFactorGrad, gradThetaGamma, preExp, ...
    gradUpsilon1, gradUpsilon2, compUpsilon1, compUpsilon2, mode, term)

% LFMGRADIENTH41 Gradient of the function h_i(z) with respect to some of the
%
%	Description:
%	hyperparameters of the kernel: m_k, C_k, D_k, m_r, C_r or D_r.
%
%	G = LFMGRADIENTH41(GAMMA1, GAMMA2, SIGMA2, GRADTHETAGAMMA, T1, T2,
%	MODE) Computes the gradient of the function h_i(z) with respect to
%	some of the parameters of the system (mass, spring or damper).
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
%	  MODE - indicates in which way the vectors t1 and t2 must be
%	   transposed
%	
%	
%
%	See also
%	LFMKERNGRADIENT, LFMXLFMKERNGRADIENT, LFMGRADIENTUPSILON


%	Copyright (c) 2007, 2008 David Luengo
%	Copyright (c) 2008 Mauricio Alvarez



% Gradient evaluation

if nargin<10,
    term = [];
end

if ~mode
    if ~term
        g = (gradUpsilon1*gradThetaGamma)*( preExp/preFactor(1) - conj(preExp)/preFactor(2)).' ...
            +(compUpsilon1*gradThetaGamma)*( - preExp/preFactorGrad(1) + conj(preExp)/preFactorGrad(2)).';
    else
        g = (gradUpsilon1*gradThetaGamma(1))*( preExp/preFactor(1)).' ...
            +(compUpsilon1*gradThetaGamma(1))*( - preExp/preFactorGrad(1)).'...
            +(conj(gradUpsilon1)*gradThetaGamma(2))*( - preExp/preFactor(2)).' ...
            +(conj(compUpsilon1)*gradThetaGamma(2))*( preExp/preFactorGrad(2)).';
    end
else
    g = (gradUpsilon1*gradThetaGamma(1))*( preExp(:,1)/preFactor(1) - preExp(:,2)/preFactor(2)).' ...
        + (compUpsilon1*gradThetaGamma(1))*( - preExp(:,1)/preFactorGrad(1) + preExp(:,2)/preFactorGrad(2)).'...
        + (gradUpsilon2*gradThetaGamma(2))*( preExp(:,2)/preFactor(4) - preExp(:,1)/preFactor(3)).' ...
        + (compUpsilon2*gradThetaGamma(2))*( - preExp(:,2)/preFactorGrad(4) + preExp(:,1)/preFactorGrad(3)).';
end
