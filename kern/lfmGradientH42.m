function    g = lfmGradientH42(preFactor, preFactorGrad, gradThetaGamma, preExp, preExpt, ...
    compUpsilon1, compUpsilon2, mode, term)

% LFMGRADIENTH42 Gradient of the function h_i(z) with respect to some of the
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
% COPYRIGHT : David Luengo, 2007, 2008
%
% COPYRIGHT : Mauricio Alvarez, 2008
%
% SEEALSO : lfmKernGradient, lfmXlfmKernGradient, lfmGradientUpsilon

% KERN


% Gradient evaluation

if nargin<9
    term = [];
end

if ~mode
    if ~term
        g = compUpsilon1*(- ( preExp/preFactorGrad(1) + preExpt/preFactor(1))*gradThetaGamma(1)...
            + ( conj(preExp)/preFactorGrad(2)+ conj(preExpt)/preFactor(2))*gradThetaGamma(2)).';
    else
        g = compUpsilon1*(- ( preExp/preFactorGrad(1) + preExpt/preFactor(1))*gradThetaGamma).'...
            + conj(compUpsilon1)*(( preExp/preFactorGrad(2)+ preExpt/preFactor(2))*gradThetaGamma).';
    end
else
    g = compUpsilon1*(( preExp(:,2)/preFactorGrad(3)+ preExpt(:,2)/preFactor(3))*gradThetaGamma(2)...
        - ( preExp(:,1)/preFactorGrad(1)+ preExpt(:,1)/preFactor(1))*gradThetaGamma(1)).'...
        - compUpsilon2*(( preExp(:,2)/preFactorGrad(4)+ preExpt(:,2)/preFactor(4))*gradThetaGamma(2)...
        - ( preExp(:,1)/preFactorGrad(2)+ preExpt(:,1)/preFactor(2))*gradThetaGamma(1)).';
end
