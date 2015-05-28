function    g = lfmGradientH(gamma1, gamma2, sigma2, gradThetaGamma, t1,t2,...
    mode, preComputeH, preComputeUpsilon)

% LFMGRADIENTH Gradient of the function h_i(z) with respect to some of the
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

if gradThetaGamma(1) == 0
    if mode==1,
        % t1 is really t1 and t2 is really t2
        Tt1 = repmat(t1, 1, size(t2, 1));
    else
        % t1 is really t2 and t2 is really t1
        %Tt1 = repmat(t2', size(t1, 1), 1);
        Tt1 = repmat(t1', size(t2, 1), 1);
    end
    g = (gradThetaGamma(2)*Tt1.*exp(-gamma2*Tt1).* preComputeUpsilon ...
        - sum(gradThetaGamma)*preComputeH)/(gamma1+gamma2);
elseif gradThetaGamma(2) == 0
    if mode==1,
       % t1 is really t1 and t2 is really t2
       Tt1 = repmat(t1, 1, size(t2, 1));
         g = (lfmGradientUpsilon(gamma1,sigma2,gradThetaGamma(1),t2, t1,2) - exp(-gamma2*Tt1) ...
        .* lfmGradientUpsilon(gamma1,sigma2,gradThetaGamma(1),t2,zeros(size(t1)),4) ...
        - sum(gradThetaGamma)*preComputeH)/(gamma1+gamma2);
    else
        % t1 is really t2 and t2 is really t1
       Tt1 = repmat(t1', size(t2, 1),1);
       %Tt1 = repmat(t2', size(t1, 1), 1);
         g = (lfmGradientUpsilon(gamma1,sigma2,gradThetaGamma(1),t2, t1,1) - exp(-gamma2*Tt1) ...
        .* lfmGradientUpsilon(gamma1,sigma2,gradThetaGamma(1),t2,zeros(size(t1)),3) ...
        - sum(gradThetaGamma)*preComputeH)/(gamma1+gamma2);
    end
end
end




