function upsilon = lfmjpComputeUpsilonVector(gamma, sigma2, t)

% LFMJPCOMPUTEUPSILONVECTOR Upsilon vector jolt. pos. with t1, t2 limits
% FORMAT
% DESC computes a portion of the LFM kernel.
% ARG gamma : Gamma value for system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% derivative wrt t1, mode 1 derivative wrt t2)
% RETURN upsilon : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio Alvarez, 2010

% KERN

sigma = sqrt(sigma2);

upsilon = gamma^2*lfmvpComputeUpsilonVector(gamma, sigma2, t) ...
    + (4/(sqrt(pi)*sigma^3))*exp(-(t.^2)./sigma2).*(t.*(gamma + 2*t/sigma2) - 1);
