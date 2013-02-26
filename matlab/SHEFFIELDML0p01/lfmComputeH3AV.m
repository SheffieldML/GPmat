function [h, compUpAV, compUpVP, compUp] =  lfmComputeH3AV(gamma1_p, gamma1_m, sigma2, t1, ...
    t2, preFactor, mode)

% LFMCOMPUTEH3AV Helper function for computing part of the LFMAV kernel.
%
%	Description:
%
%	H = LFMCOMPUTEH3AV(GAMMA1, GAMMA2, SIGMA2, T1, T2, PREFACTOR, MODE)
%	computes a portion of the LFMAV kernel.
%	 Returns:
%	  H - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA1 - Gamma value for first system.
%	  GAMMA2 - Gamma value for second system.
%	  SIGMA2 - length scale of latent process.
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1).
%	  PREFACTOR - precomputed constants.
%	  MODE - indicates the correct derivative.


%	Copyright (c) 2010 Mauricio Alvarez


% Evaluation of h

if nargout>1
    [compUpAV{1}, compUpVP{1}, compUp{1}] = lfmavComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode);
    [compUpAV{2}, compUpVP{2}, compUp{2}] = lfmavComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
    h = preFactor(1)*compUpAV{1} + preFactor(2)*compUpAV{2};
else
    h = preFactor(1)*lfmavComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode) ...
        + preFactor(2)*lfmavComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
end
