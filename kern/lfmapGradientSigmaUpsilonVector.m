function dUpsilonS = lfmapGradientSigmaUpsilonVector(gamma, sigma2, t)

% LFMAPGRADIENTSIGMAUPSILONVECTOR Gradient of upsilon vector ap wrt sigma
% FORMAT
% DESC computes the gradient of a portion of the LFMAP kernel.
% ARG gamma : Gamma value for system.
% ARG sigma2 : length scale of latent process.
% ARG t : first time input (number of time points x 1).
% RETURN upsilon : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio Alvarez, 2010
%
% SEEALSO : lfmapComputeUpsilonMatrix.m

% KERN


dUpsilon = lfmGradientSigmaUpsilonVector(gamma, sigma2, t);

dUpsilonS = (gamma^2)*dUpsilon + (2/(sqrt(pi)*sigma2))* ...
    exp(-(t.^2)./sigma2).*(gamma + 2*t/sigma2).*(1-2*(t.^2)/sigma2) ...
    + (8/(sqrt(pi)*sigma2^2))*t.*exp(-(t.^2)./sigma2);
