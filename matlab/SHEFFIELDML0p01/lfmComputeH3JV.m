function [h, compUpJV] =  lfmComputeH3JV(gamma1_p, gamma1_m, sigma2, t1, ...
    t2, preFactor, mode)

% LFMCOMPUTEH3JV Helper function for computing part of the LFMJV kernel.
%
%	Description:
%
%	H = LFMCOMPUTEH3JV(GAMMA1, GAMMA2, SIGMA2, T1, T2, PREFACTOR, MODE)
%	computes a portion of the LFMJV kernel.
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
    compUpJV{1} = lfmjvComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode);
    compUpJV{2} = lfmjvComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
    h = preFactor(1)*compUpJV{1} + preFactor(2)*compUpJV{2};
else
    h = preFactor(1)*lfmjvComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode) ...
        + preFactor(2)*lfmjvComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
end
