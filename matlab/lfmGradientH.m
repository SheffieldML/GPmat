function    g = lfmGradientH(gamma1, gamma2, sigma2, gradThetaGamma, t1, t2)

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
% ARG t2 : second time input (number of time points x 1).
% RETURN g : Gradient of the function with respect to the desired
% parameter.
%
% COPYRIGHT : David Luengo, 2007
%
% MODIFICATIONS : David Luengo, 2008
%
% SEEALSO : lfmKernGradient, lfmXlfmKernGradient, lfmGradientUpsilon

% LFM


% Creation of the time matrices

Tt1 = repmat(t1,1,size(t2, 1));
Tt2 = repmat(t2',size(t1, 1),1);

% Gradient evaluation

g = (lfmGradientUpsilon(gamma1,sigma2,gradThetaGamma(1),Tt2,Tt1) ...
    + gradThetaGamma(2)*Tt1.*exp(-gamma2*Tt1) ...
    .* lfmComputeUpsilon(gamma1,sigma2,Tt2,zeros(size(Tt1))) - exp(-gamma2*Tt1) ...
    .* lfmGradientUpsilon(gamma1,sigma2,gradThetaGamma(1),Tt2,zeros(size(Tt1))) ...
    - sum(gradThetaGamma)*lfmComputeH(gamma1,gamma2,sigma2,t1,t2))/(gamma1+gamma2);
