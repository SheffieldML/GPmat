function    g = lfmGradientH42VP(preFactor, preFactorGrad, gradThetaGamma, ...
    preExp, preExpg, preExpt, compUpsilon1, compUpsilon2)

% LFMGRADIENTH42VP Gradient of the function h_i(z) with respect to some of the
% hyperparameters of the kernel: m_k, C_k, D_k, m_r, C_r or D_r.
% FORMAT
% DESC Computes the gradient of the function h_i(z) with respect to some of
% the parameters of the system (mass, spring or damper).
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG sigma2 : length scale of latent process.
% ARG gradThetaGamma : Vector with the gradient of gamma1 and gamma2 with
% respect to the desired parameter.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1)
% RETURN g : Gradient of the function with respect to the desired
% parameter.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010
%
% SEEALSO : lfmKernGradient, lfmXlfmKernGradient, lfmGradientUpsilon

% KERN


% Gradient evaluation


g = compUpsilon1*(( - preExpg(:,2)/preFactorGrad(3) - preExpt(:,2)/preFactor(3) + preExp(:,2)/preFactor(3))*gradThetaGamma(2)...
    + ( preExpg(:,1)/preFactorGrad(1) + preExpt(:,1)/preFactor(1) - preExp(:,1)/preFactor(1))*gradThetaGamma(1)).'...
    - compUpsilon2*(( - preExpg(:,2)/preFactorGrad(4) - preExpt(:,2)/preFactor(4) + preExp(:,2)/preFactor(4))*gradThetaGamma(2)...
    + ( preExpg(:,1)/preFactorGrad(2) + preExpt(:,1)/preFactor(2) - preExp(:,1)/preFactor(2))*gradThetaGamma(1)).';

    

