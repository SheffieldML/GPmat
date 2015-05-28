function [h, compUpVV, compUpVP, compUp] =  lfmComputeH3VV(gamma1_p, gamma1_m, sigma2, t1, ...
    t2, preFactor, mode)

% LFMCOMPUTEH3VV Helper function for computing part of the LFMVXLFMV kernel.
% FORMAT
% DESC computes a portion of the LFMVXLFMV kernel.
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG preFactor : precomputed constants.
% ARG mode: indicates the correct precomputations.
% RETURN h : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

% Evaluation of h

if nargout>1
    [compUpVV{1}, compUpVP{1}, compUp{1}] = lfmvvComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode);
    [compUpVV{2}, compUpVP{2}, compUp{2}] = lfmvvComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
    h = preFactor(1)*compUpVV{1} + preFactor(2)*compUpVV{2};
else
    h = preFactor(1)*lfmvvComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode) ...
        + preFactor(2)*lfmvvComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);  
    
%     %h = preFactor(2)*lfmvvComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);  
%     
%      h = lfmvvComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
% 
%      epsilon = 1e-6;
%      valg = gamma1_m;
%      gamma1_m = valg + epsilon;
%      h1 = lfmvvComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
%      gamma1_m = valg - epsilon;
%      h2 = lfmvvComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
%      gamma1_m = valg;
%      numerics = 0.5*(h1-h2)/epsilon;
%      theo = lfmvvGradientUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
     
     
end
