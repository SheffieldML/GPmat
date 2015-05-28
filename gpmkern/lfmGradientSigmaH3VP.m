function    g = lfmGradientSigmaH3VP(gamma1, gamma2, sigma2, t1, t2, ...
    preFactor, mode)

% LFMGRADIENTSIGMAH3VP Gradient of the function h_i(z) with respect \sigma.
% FORMAT
% DESC Computes the gradient of the function h_i(z) with respect to the
% length-scale of the input "force", \sigma.
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG mode: indicates in which way the vectors t1 and t2 must be transposed
% RETURN g : Gradient of the function with respect to \sigma.
%
% COPYRIGHT : Mauricio Alvarez, 2010
%
% SEEALSO : lfmKernGradient, lfmXlfmKernGradient, lfmGradientSigmaUpsilon

% KERN

g =  preFactor(1)*lfmvpGradientSigmaUpsilonMatrix(gamma1,sigma2,t1,t2, mode) + ...
    preFactor(2)*lfmvpGradientSigmaUpsilonMatrix(gamma2,sigma2,t1,t2, mode);


