function upsilon = lfmvpComputeUpsilonDiagVector(gamma, sigma2, t, mode)

% LFMVPCOMPUTEUPSILONDIAGVECTOR Upsilon diag vector vel. pos. with t1, t2 limits
% FORMAT
% DESC computes a portion of the LFMVP kernel.
% ARG gamma : Gamma value for system.
% ARG sigma2 : length scale of latent process.
% ARG t : first time input (number of time points x 1).
% ARG mode : operation mode, according to the derivative (mode 0,
% derivative wrt t1, mode 1 derivative wrt t2)
% RETURN upsilon : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio Alvarez, 2010
%
% SEEALSO : lfmComputeUpsilonMatrix.F, lfmComputeH3.m

% KERN

sigma = sqrt(sigma2);

if mode==0   
    upsilon = -gamma*lfmComputeUpsilonDiagVector(gamma, sigma2, t) ...
        + (2/(sqrt(pi)*sigma));    
else
    upsilon = gamma*lfmComputeUpsilonDiagVector(gamma, sigma2, t) ...
        - (2/(sqrt(pi)*sigma)) ...
        + (2/(sqrt(pi)*sigma))*exp(-gamma*t).*(exp(-(t.^2)/sigma2));
end
