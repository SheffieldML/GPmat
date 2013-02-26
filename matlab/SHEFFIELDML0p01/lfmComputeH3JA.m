function [h, compUp] =  lfmComputeH3JA(gamma1_p, gamma1_m, sigma2, t1, ...
    t2, preFactor, mode)

% LFMCOMPUTEH3JA Helper function for computing part of the LFMJA kernel.
%
%	Description:
%
%	H = LFMCOMPUTEH3JA(GAMMA1, GAMMA2, SIGMA2, T1, T2, PREFACTOR, MODE)
%	computes a portion of the LFMJA kernel.
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
    compUp = cell(2,1);
    compUp{1} = lfmjaComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode);
    compUp{2} = lfmjaComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
    h = preFactor(1)*compUp{1} + preFactor(2)*compUp{2};
else
    h = preFactor(1)*lfmjaComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode) ...
        + preFactor(2)*lfmjaComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
end
