function [h, compUpJP] =  lfmComputeH3JP(gamma1_p, gamma1_m, sigma2, t1, ...
    t2, preFactor, mode)

% LFMCOMPUTEH3JP Helper function for computing part of the LFMJP kernel.
% FORMAT
% DESC computes a portion of the LFMAP kernel.
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG preFactor : precomputed constants.
% ARG mode: indicates the correct derivative.
% RETURN h : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio Alvarez, 2010

% KERN

% Evaluation of h

if nargout>1    
    compUpJP{1} = lfmjpComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode);
    compUpJP{2} = lfmjpComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
    h = preFactor(1)*compUpJP{1} + preFactor(2)*compUpJP{2};
else
    h = preFactor(1)*lfmjpComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode) ...
        + preFactor(2)*lfmjpComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
end
