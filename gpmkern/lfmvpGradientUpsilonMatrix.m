function dUpsilon = lfmvpGradientUpsilonMatrix(gamma, sigma2, t1, ...
    t2, mode, upsilon)

% LFMVPGRADIENTUPSILONMATRIX Gradient upsilon matrix vel. pos.
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

sigma = sqrt(sigma2);

if nargin<6
    upsilon = lfmComputeUpsilonMatrix(gamma, sigma2, t1, t2);
    if nargin <5
        mode =0;
    end
end

if mode ==0
    dUpsilon = -upsilon - gamma*lfmGradientUpsilonMatrix(gamma, sigma2, t1, t2);
else
    dUpsilon = upsilon + gamma*lfmGradientUpsilonMatrix(gamma, sigma2, t1, t2) ...
        - (2/(sqrt(pi)*sigma))*(t1.*exp(-gamma*t1))*(exp(-(t2.^2)/sigma2)).';
end

