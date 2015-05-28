function dUpsilon = lfmvpGradientUpsilonVector(gamma, sigma2, t, upsilon)

% LFMVPGRADIENTUPSILONVECTOR Gradient upsilon vector vel. pos.
% FORMAT
% DESC computes the gradient of a portion of the LFM kernel.
% ARG gamma : Gamma value for system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG upsilon : precomputation of the upsilon matrix.
% ARG mode : operation mode, according to the derivative (mode 0,
% derivative wrt t1, mode 1 derivative wrt t2)
% RETURN upsilon : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010
%
% SEEALSO : lfmvpComputeUpsilonMatrix.m

% KERN

if nargin<4
    upsilon = lfmComputeUpsilonVector(gamma, sigma2, t);
end

dUpsilon = -upsilon - gamma*lfmGradientUpsilonVector(gamma, sigma2, t);

