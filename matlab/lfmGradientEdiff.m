function    matGrad = lfmGradientEdiff(gamma,sigma2,gradThetaGamma,Tt,Tt2);
%
%
% Gradient of E_r(gamma_k,t,t') with respect to theta_k, where theta_k can
% be m_k, C_k or D_k.
%
%   matGrad = GPLfmGradThetaE(gamma,sigma2,gradThetaGamma,Tt,Tt2);
%
%   matGrad        - Matrix with the gradients of the kernel.
%   gamma          - Propagation constant of the system.
%   sigma2         - Scalar. Length scale of the r-th input.
%   gradThetaGamma - Gradient of gamma with respect to theta_k.
%   Tt             - Matrix TxT2. Time instants for x_k.
%   Tt2            - Matrix TxT2. Time instants for f_r.
%
% Author            : David Luengo Garcia
% Place and Date    : Manchester, 26 October 2007
% Last Modification : 29 October 2007


% Parameters of the function

sigma = sqrt(sigma2);
% Gradient
matGrad = (sigma*gradThetaGamma/2)*...
    (exp(-((Tt2+sigma2*gamma/2).^2)/sigma2) - ...
    exp(-((Tt-(Tt2+sigma2*gamma/2)).^2)/sigma2));

