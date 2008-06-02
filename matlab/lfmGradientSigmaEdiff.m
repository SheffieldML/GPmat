function    matGrad = lfmGradientSigmaEdiff(gamma,sigma2,Tt,Tt2);
%
%
% Gradient of E_r(gamma,T,T') with respect to sigma_r.
%
%   matGrad = GPLfmGradSigmaE(gamma,sigma2,Tt,Tt2);
%
%   matGrad - Matrix with the gradients of the function E_r.
%   gamma   - Propagation constant of the system.
%   sigma2  - Scalar. Length scale of the r-th input.
%   Tt      - Matrix TxT2. Time instants for x_k.
%   Tt2     - Matrix TxT2. Time instants for f_r.
%
% Author            : David Luengo Garcia
% Place and Date    : Manchester, 26 October 2007
% Last Modification : 26 October 2007


% Parameters of the function

% Gradient

matGrad = -(((Tt-Tt2)/sigma2+gamma/2).* ...
    exp(-((Tt-(Tt2+sigma2*gamma/2)).^2)/sigma2) + (Tt2/sigma2-gamma/2) ...
    .*exp(-((Tt2+sigma2*gamma/2).^2)/sigma2));

