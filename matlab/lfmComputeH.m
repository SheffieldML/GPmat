function  h = lfmComputeH(gamma1, gamma2, sigma2, t1, t2);

% LFMCOMPUTEH Helper function for comptuing part of the LFM kernel.
% FORMAT
% DESC computes a portion of the LFM kernel.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG l : length scale of latent process.
% RETURN h : result of this subcomponent of the kernel for the given values.
%
% DESC computes a portion of the LFM kernel and gradients with
% respect to various parameters.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG l : length scale of latent process.
% RETURN h : result of this subcomponent of the kernel for the given values.
% RETURN grad_D_gamma1 : gradient of H with respect to GAMMA1.
% RETURN grad_D_gamma2 : gradient of H with respect to GAMMA1.
% RETURN grad_L : gradient of H with respect to length scale of
% latent process.
%
% COPYRIGHT : David Luengo, 2007
%
% MODIFICATIONS : Neil D. Lawrence, 2006
%
% SEEALSO : lfmKernParamInit

% KERN

sigma = sqrt(sigma2);

% Creation of the time matrices
Tt1 = repmat(t1, 1, size(t2, 1));
Tt2 = repmat(t2', size(t1, 1), 1);

% Evaluation of h
h = (exp(sigma2*(gamma1^2)/4)/(gamma1+gamma2))* ...
     (exp(gamma1*(Tt1-Tt2)).*lfmComputeEdiff(gamma1,sigma2,Tt2,Tt1) ...
      - exp(-(gamma1*Tt2+gamma2*Tt1)).*lfmComputeEdiff(gamma1,sigma2,Tt2,0));



function x = lfmComputeEdiff(gamma,sigma2,Tt,Tt2);
%
%
% Function that evaluates the function E_i, used in Kxf and Kxx. Requires
% the function erfz, which provides the erf(z) for complex values.
%
%   Hi = GPLfmEvalEi(gamma,sigma2,Tt,Tt2);
%
%   Ei     - Matrix TxT2 with the desired values of E_i.
%   gamma  - Propagation constant.
%   sigma2 - Length scale of the GP input.
%   Tt     - Matrix TxT2. Time instants for x_k.
%   Tt2    - Matrix TxT2. Time instants for x_r.
%
% Author            : David Luengo Garcia
% Place and Date    : Manchester, 29 October 2007
% Last Modification : 29 October 2007


% Parameters of the kernel

sigma = sqrt(sigma2);

% Initialization of vectors and matrices

x = zeros(size(Tt));

% Evaluation of Ei

x = erfz((Tt-Tt2)/sigma-sigma*gamma/2) + erfz(Tt2/sigma+sigma*gamma/2);

return;
