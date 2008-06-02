function    matGrad = lfmGradientH(gamma1, gamma2, sigma2, gradThetaGamma, t, t2)

%
% Gradient of h_i(gamma1,gamma2,t,t') w.r.t. theta_l, where theta_l can be
% m_k, C_k, D_k, m_r, C_r or D_r.
%
%   matGrad = GPLfmgradThetaHi(gamma1,gamma2,sigma2,gradThetaGamma,t,t2);
%
%   matGrad        - Matrix with the gradients of H_i.
%   gamma1         - First value of gamma (gamma_r).
%   gamma2         - Second value of gamma (gamma_k).
%   sigma2         - Vector Mx1. Length scale of the M input "forces".
%   gradThetaGamma - Vector 2x1. Gradient of gamma1 and gamma2 with respect
%                    to the desired parameter.
%   t              - Vector Tx1. Time instants for x_k.
%   t2             - Vector T2x1. Time instants for x_r.
%
% Author            : David Luengo Garcia
% Place and Date    : Manchester, 28 October 2007
% Last Modification : 7 November 2007


% Creation of the time matrices
Tt = repmat(t,1,size(t2, 1));
Tt2 = repmat(t2',size(t, 1),1);

% Gradient evaluation

gradThetaGammai = gradThetaGamma(1)*(Tt-Tt2) ...
    .*exp(gamma1*(Tt-Tt2)).*lfmComputeEdiff(gamma1,sigma2,Tt2,Tt) ...
    + exp(gamma1*(Tt-Tt2)).* ...
    lfmGradientEdiff(gamma1,sigma2,gradThetaGamma(1),Tt2,Tt) ...
    + (gradThetaGamma(1)*Tt2 + gradThetaGamma(2)*Tt) ...
    .*exp(-(gamma1*Tt2+gamma2*Tt)).*lfmComputeEdiff(gamma1,sigma2,Tt2,0) ...
    - exp(-(gamma1*Tt2+gamma2*Tt)).* ...
    lfmGradientEdiff(gamma1,sigma2,gradThetaGamma(1),Tt2,0);

matGrad = ((sigma2*gamma1/2)*gradThetaGamma(1) - ...
           (sum(gradThetaGamma)/(gamma1+gamma2))) * ...
          lfmComputeH(gamma1,gamma2,sigma2,t,t2) + ...
          (exp(sigma2*(gamma1^2)/4)/(gamma1+gamma2)) * ...
          gradThetaGammai;

