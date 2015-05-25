function dUpsilon = lfmvvGradientUpsilonMatrix(gamma, sigma2, t1, ...
    t2, mode, upsilon)

% LFMVVGRADIENTUPSILONMATRIX Gradient upsilon matrix vel. vel.
% FORMAT
% DESC computes the gradient of a portion of the LFMVV kernel.
% ARG gamma : Gamma value for system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG upsilon : precomputation of the upsilon matrix.
% ARG mode : operation mode, according to the derivative (mode 0,
% derivative wrt t1, mode 1 derivative wrt t2)
% RETURN upsilon : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio Alvarez, 2010
%
% SEEALSO : lfmvvComputeUpsilonMatrix.m

% KERN

sigma = sqrt(sigma2);

if nargin<6
    upsilon = lfmvpComputeUpsilonMatrix(gamma, sigma2, t1, t2, mode);
    if nargin <5
        mode =0;
    end
end

if mode ==0
    dUpsilon = upsilon + gamma*lfmvpGradientUpsilonMatrix(gamma, sigma2, t1, t2, mode) ...
        - (2/(sqrt(pi)*sigma))*((1-gamma*t1).*exp(-gamma*t1))*(exp(-(t2.^2)/sigma2)).';
else
    dUpsilon = -upsilon - gamma*lfmvpGradientUpsilonMatrix(gamma, sigma2, t1, t2, mode);
end

