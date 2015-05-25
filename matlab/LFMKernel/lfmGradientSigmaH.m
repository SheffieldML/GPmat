function    g = lfmGradientSigmaH(gamma1, gamma2, sigma2, t1, t2);

% LFMGRADIENTSIGMAH Gradient of the function h_i(z) with respect \sigma.
% FORMAT
% DESC Computes the gradient of the function h_i(z) with respect to the
% length-scale of the input "force", \sigma.
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% RETURN g : Gradient of the function with respect to \sigma.
%
% COPYRIGHT : David Luengo, 2007
%
% MODIFICATIONS : David Luengo, 2008
%
% SEEALSO : lfmKernGradient, lfmXlfmKernGradient, lfmGradientSigmaUpsilon

% LFM


% Creation of the time matrices

Tt1 = repmat(t1,1,size(t2, 1));
Tt2 = repmat(t2',size(t1, 1),1);

% Gradient

g = (lfmGradientSigmaUpsilon(gamma1,sigma2,Tt2,Tt1) - exp(-gamma2*Tt1) ...
    .*lfmGradientSigmaUpsilon(gamma1,sigma2,Tt2,zeros(size(Tt1))))/(gamma1+gamma2);
