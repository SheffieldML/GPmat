function g = lfmGradientH41VP(preFactor, preFactorGrad, gradThetaGamma, preExp, gradUpsilon1, gradUpsilon2, compUpsilon1, compUpsilon2)

% LFMGRADIENTH41VP Gradient of the function h_i(z) with respect to some of the
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
% ARG mode: indicates in which way the vectors t1 and t2 must be transposed
% RETURN g : Gradient of the function with respect to the desired
% parameter.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN


% Gradient evaluation

g = (gradUpsilon1*gradThetaGamma(1))*(preExp(:,2)/preFactor(2) - preExp(:,1)/preFactor(1)).' ...
    + (compUpsilon1*gradThetaGamma(1))*( preExp(:,1)/preFactorGrad(1) - preExp(:,2)/preFactorGrad(2)).'...
    + (gradUpsilon2*gradThetaGamma(2))*( preExp(:,1)/preFactor(3) - preExp(:,2)/preFactor(4)).' ...
    + (compUpsilon2*gradThetaGamma(2))*( preExp(:,2)/preFactorGrad(4) - preExp(:,1)/preFactorGrad(3)).';

