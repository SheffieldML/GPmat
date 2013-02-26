function  h = lfmComputeH(gamma1, gamma2, sigma2, t1, t2, mode)

% LFMCOMPUTEH Helper function for computing part of the LFM kernel.
%
%	Description:
%
%	H = LFMCOMPUTEH(GAMMA1, GAMMA2, SIGMA2, T1, T2, MODE) computes a
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
%	
%
%	See also
%	LFMKERNPARAMINIT, LFMXLFMKERNCOMPUTE


%	Copyright (c) 2007, Mauricio Alvarez, 2008 David Luengo


%	With modifications by Neil D. Lawrence 2007


% Evaluation of h


if mode==1
    % t1 is really t1 and t2 is really t2
%    Tt1 = repmat(t2', size(t1, 1), 1);
    Tt1 = repmat(t1, 1, size(t2, 1));
    h = (lfmComputeUpsilon(gamma1,sigma2,t2,t1,1) - exp(-gamma2*Tt1) ...
        .* lfmComputeUpsilon(gamma1,sigma2,t2, zeros(size(t1)),3))/(gamma1+gamma2);
else
    % t1 is really t2 and t2 is really t1
    Tt1 = repmat(t1', size(t2, 1), 1);
    %Tt1 = repmat(t2', size(t1, 1), 1);
%    Tt1 = repmat(t1, 1, size(t2, 1))';
    h = (lfmComputeUpsilon(gamma1,sigma2,t2,t1,2) - exp(-gamma2*Tt1) ...
        .* lfmComputeUpsilon(gamma1,sigma2,t2, zeros(size(t1)),4))/(gamma1+gamma2);
end

