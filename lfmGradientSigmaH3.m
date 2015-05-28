function    g = lfmGradientSigmaH3(gamma1, gamma2, sigma2, t1, t2, preFactor, mode, term)

% LFMGRADIENTSIGMAH3 Gradient of the function h_i(z) with respect \sigma.
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

if nargin<8
    term = [];
end

if ~mode
    if ~term
        g = preFactor*lfmGradientSigmaUpsilonMatrix(gamma1,sigma2, t1,t2);
    else
        gradupsilon = lfmGradientSigmaUpsilonMatrix(gamma1,sigma2, t1,t2);
        g = -preFactor(1)*gradupsilon + preFactor(2)*conj(gradupsilon);
    end
else
    g =  preFactor(1)*lfmGradientSigmaUpsilonMatrix(gamma1,sigma2,t1,t2) + ...
        preFactor(2)*lfmGradientSigmaUpsilonMatrix(gamma2,sigma2,t1,t2);
end



