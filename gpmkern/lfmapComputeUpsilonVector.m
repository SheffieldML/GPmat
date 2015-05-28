function [upsilonap, upsilon] = lfmapComputeUpsilonVector(gamma, sigma2, t1)

% LFMAPCOMPUTEUPSILONVECTOR Upsilon vector for acce. pos. with t1 limit
% FORMAT
% DESC computes a portion of the LFMAP kernel.
% ARG gamma : Gamma value for system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% RETURN upsilon : result of this subcomponent of the kernel for the given
% values.
%
% COPYRIGHT : Mauricio Alvarez, 2010
%
% SEEALSO : lfmapComputeUpsilonMatrix.m

% KERN

sigma = sqrt(sigma2);

if nargout > 1
    upsilon = lfmComputeUpsilonVector(gamma, sigma2, t1);
    upsilonap = gamma^2*upsilon - (2/(sqrt(pi)*sigma))*exp(-(t1.^2)/sigma2).*(gamma + 2*t1/sigma2);
else
    upsilonap = gamma^2*lfmComputeUpsilonVector(gamma, sigma2, t1) ...
        - (2/(sqrt(pi)*sigma))*exp(-(t1.^2)/sigma2).*(gamma + 2*t1/sigma2);
end
