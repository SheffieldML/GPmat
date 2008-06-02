function    matGrad = lfmGradientUpsilon(gamma,sigma2,gradThetaGamma,t,t2);
%
%
% Gradient of upsilon_r(gamma_k,t,t') with respect to theta_k, where theta_k
% can be m_k, C_k or D_k.
%
%   matGrad = GPLfmgradThetaupsilon(gamma,sigma2,gradThetaGamma,t,t2);
%
%   matGrad        - Matrix with the gradients of the kernel.
%   gamma          - Propagation constant of the system.
%   sigma2         - Scalar. Length scale of the r-th input.
%   gradThetaGamma - Gradient of gamma with respect to theta_k.
%   t              - Vector Tx1. Time instants for x_k.
%   t2             - Vector T2x1. Time instants for f_r.
%
% Author            : David Luengo Garcia
% Place and Date    : Manchester, 26 October 2007
% Last Modification : 26 October 2007


%Parameters of the function

T = length(t);
T2 = length(t2);

sigma = sqrt(sigma2);

% Initialization of vectors and matrices

Tt = zeros(T,T2);
Tt2 = zeros(T,T2);

upsilon = zeros(T,T2);
matGrad = zeros(T,T2);

% Creation of the time matrices

t = reshape(t,T,1);
Tt = repmat(t,1,T2);
t2 = reshape(t2,1,T2);
Tt2 = repmat(t2,T,1);

% Gradient

upsilon = exp(sigma2*(gamma^2)/4)*exp(-gamma*(Tt-Tt2))...
    .*(erfz((Tt-Tt2)/sigma-sigma*gamma/2) + erfz(Tt2/sigma+sigma*gamma/2));

matGrad = gradThetaGamma*(sigma2*gamma/2-(Tt-Tt2)).*upsilon + ...
    exp(sigma2*(gamma^2)/4)*exp(-gamma*(Tt-Tt2)).*...
    GPLfmgradThetaE(gamma,sigma2,gradThetaGamma,Tt,Tt2);

return;
