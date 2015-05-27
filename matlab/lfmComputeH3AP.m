function [h, compUpAP, compUp] =  lfmComputeH3AP(gamma1_p, gamma1_m, sigma2, t1, ...
    t2, preFactor, mode)

% LFMCOMPUTEH3AP Helper function for computing part of the LFMAP kernel.
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
    [compUpAP{1}, compUp{1}] = lfmapComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode);
    [compUpAP{2}, compUp{2}] = lfmapComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
    h = preFactor(1)*compUpAP{1} + preFactor(2)*compUpAP{2};
else
    h = preFactor(1)*lfmapComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode) ...
        + preFactor(2)*lfmapComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
end
