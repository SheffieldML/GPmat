function    g = lfmGradientSigmaH(gamma1, gamma2, sigma2, t1, t2, mode)

% LFMGRADIENTSIGMAH Gradient of the function h_i(z) with respect \sigma.
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
% COPYRIGHT : David Luengo, 2007, 2008, Mauricio Alvarez, 2008
%
% MODIFICATIONS : Mauricio Alvarez, 2008

% SEEALSO : lfmKernGradient, lfmXlfmKernGradient, lfmGradientSigmaUpsilon

% KERN


% Gradient


if mode==1,
    % t1 is really t1 and t2 is really t2
    Tt1 = repmat(t1, 1, size(t2, 1));
    g = (lfmGradientSigmaUpsilon(gamma1,sigma2,t2,t1,2) - exp(-gamma2*Tt1) ...
        .*lfmGradientSigmaUpsilon(gamma1,sigma2,t2,zeros(size(t1,1)),4))/(gamma1+gamma2);
else
    % t1 is really t2 and t2 is really t1
    Tt1 = repmat(t2', size(t1, 1), 1);
    g = (lfmGradientSigmaUpsilon(gamma1,sigma2,t2,t1,1) - exp(-gamma2*Tt1) ...
        .*lfmGradientSigmaUpsilon(gamma1,sigma2,t2,zeros(size(t1,1)),3))/(gamma1+gamma2);
end

