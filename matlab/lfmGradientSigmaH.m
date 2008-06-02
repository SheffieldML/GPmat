function    MatGrad = lfmGradientSigmaH(gamma1,gamma2,sigma2,t,t2);
%
%
% Gradient of H_i(gamma1,gamma2,t,t') with respect to sigma_i.
%
%   MatGrad = GPLfmGradSigmaHi(gamma1,gamma2,sigma2,t,t2);
%
%   MatGrad - Matrix with the gradients of the function H_i.
%   gamma1  - First value of gamma (gamma_r).
%   gamma2  - Second value of gamma (gamma_k).
%   sigma2  - Scalar. Length scale of the i-th input.
%   t       - Vector Tx1. Time instants for x_k.
%   t2      - Vector T2x1. Time instants for x_r.
%
% Author            : David Luengo Garcia
% Place and Date    : Manchester, 27 October 2007
% Last Modification : 28 October 2007


%Parameters of the function


sigma = sqrt(sigma2);



% Creation of the time matrices

Tt = repmat(t, 1, size(t2, 1));
Tt2 = repmat(t2', size(t, 1), 1);

% Gradient

Hi = lfmComputeH(gamma1,gamma2,sigma2,t,t2);

MatGrad = (sigma*(gamma1^2)/2)*Hi + (exp(sigma2*(gamma1^2)/4)/(gamma1+gamma2)) ...
    *(exp(gamma1*(Tt-Tt2)).*lfmGradientSigmaEdiff(gamma1,sigma2,Tt2,Tt) ...
    - exp(-(gamma1*Tt2+gamma2*Tt)).*lfmGradientSigmaEdiff(gamma1,sigma2,Tt2,0));


