function dUpsilon = lfmapGradientUpsilonVector(gamma, sigma2, t, upsilon)

% LFMAPGRADIENTUPSILONVECTOR Gradient upsilon vector accel. pos.
% FORMAT
% DESC computes the gradient of a portion of the LFM kernel.
% ARG gamma : Gamma value for system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG upsilon : precomputation of the upsilon matrix.
% RETURN upsilon : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio Alvarez, 2010
%
% SEEALSO : lfmapComputeUpsilonVector.m

% KERN

sigma = sqrt(sigma2);

if nargin<4
    upsilon = lfmComputeUpsilonVector(gamma, sigma2, t);
end

dUpsilon = 2*gamma*upsilon + gamma^2*lfmGradientUpsilonVector(gamma, sigma2, t) ...
        - (2/(sqrt(pi)*sigma))*exp(-(t.^2)./sigma2);

