function [h, compUp] =  lfmComputeH3(gamma1_p, gamma1_m, sigma2, t1,t2,preFactor,...
    mode, term)

% LFMCOMPUTEH3 Helper function for computing part of the LFM kernel.
%
%	Description:
%
%	H = LFMCOMPUTEH3(GAMMA1, GAMMA2, SIGMA2, T1, T2, MODE) computes a
%	portion of the LFM kernel.
%	 Returns:
%	  H - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA1 - Gamma value for first system.
%	  GAMMA2 - Gamma value for second system.
%	  SIGMA2 - length scale of latent process.
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1).
%	  MODE - indicates in which way the vectors t1 and t2 must be
%	   transposed
%	
%	
%
%	See also
%	LFMKERNPARAMINIT, LFMXLFMKERNCOMPUTE


%	Copyright (c) 2008 Mauricio Alvarez


%	With modifications by Neil D. Lawrence 2007


% Evaluation of h

if nargin<8
    term =[];
end

if ~mode
    if ~term
        if nargout >1
            compUp = lfmComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2);
            h = preFactor*compUp;
        else
            h = preFactor*lfmComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2);
        end
    else
        if nargout > 1
            compUp = lfmComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2);
            h = -preFactor(1)*compUp + preFactor(2)*conj(compUp);
        else
            upsilon = lfmComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2);
            h = -preFactor(1)*upsilon + preFactor(2)*conj(upsilon);
        end
    end
else
    if nargout>1
        compUp = cell(2,1);
        compUp{1} = lfmComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2);
        compUp{2} = lfmComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2);
        h = preFactor(1)*compUp{1} + preFactor(2)*compUp{2};
    else
        h = preFactor(1)*lfmComputeUpsilonMatrix(gamma1_p,sigma2, t1,t2) ...
            + preFactor(2)*lfmComputeUpsilonMatrix(gamma1_m,sigma2, t1,t2);
    end
end