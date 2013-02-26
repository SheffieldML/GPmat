function [h, compUpV, compUp] =  lfmComputeH3VP(gamma1_p, gamma1_m, sigma2, t1, ...
    t2, preFactor, mode)

% LFMCOMPUTEH3VP Helper function for computing part of the LFMVXLFM kernel.
%
%	Description:
%
%	H = LFMCOMPUTEH3VP(GAMMA1, GAMMA2, SIGMA2, T1, T2, PREFACTOR, MODE)
%	computes a portion of the LFMVXLFM kernel.
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
%	  MODE - indicates the correct precomputations.


%	Copyright (c) 2010 Mauricio Alvarez


% Evaluation of h

if nargout>1
    compUpV = cell(2,1);
    compUp = cell(2,1);
    [compUpV{1}, compUp{1}] = lfmvpComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode);
    [compUpV{2}, compUp{2}] = lfmvpComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
    h = preFactor(1)*compUpV{1} + preFactor(2)*compUpV{2};
else
    h = preFactor(1)*lfmvpComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2, mode) ...
        + preFactor(2)*lfmvpComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2, mode);
end
