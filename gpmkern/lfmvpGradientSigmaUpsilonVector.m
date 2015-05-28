function dUpsilonS = lfmvpGradientSigmaUpsilonVector(gamma, sigma2, t)

% LFMVPGRADIENTSIGMAUPSILONVECTOR Gradient of upsilon vector vp wrt sigma
% FORMAT
% DESC computes the gradient of a portion of the LFMVP kernel.
% ARG gamma : Gamma value for system.
% ARG sigma2 : length scale of latent process.
% ARG t : first time input (number of time points x 1).
% RETURN upsilon : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio Alvarez, 2010
%
% SEEALSO : lfmvpComputeUpsilonMatrix.m

% KERN


dUpsilon = lfmGradientSigmaUpsilonVector(gamma, sigma2, t);

dUpsilonS = -gamma*dUpsilon-(2/(sqrt(pi)*sigma2))*(1-2*(t.^2)/sigma2).* ...
    exp(-(t.^2)./sigma2);
