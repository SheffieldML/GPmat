function  h = lfmwhiteComputeH(gamma1, gamma2, t1, t2, isStationary1, isStationary2)

% LFMWHITECOMPUTEH Helper function for computing part of the LFM-WHITE
%
%	Description:
%	kernel.
%
%	H = LFMWHITECOMPUTEH(GAMMA1, GAMMA2, T1, T2, ISSTATIONARY1,
%	ISSTATIONARY2) computes a portion of the LFM-WHITE kernel.
%	 Returns:
%	  H - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA1 - Gamma value for first system.
%	  GAMMA2 - Gamma value for second system.
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1).
%	  ISSTATIONARY1 - indicates whether the stationary version of the
%	   first kernel is used (1) or not (0).
%	  ISSTATIONARY2 - indicates whether the stationary version of the
%	   second kernel is used (1) or not (0).
%	
%
%	See also
%	LFMWHITEKERNPARAMINIT, LFMWHITEXLFMWHITEKERNCOMPUTE


%	Copyright (c) 2009 David Luengo



if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1 - T2;
indT = deltaT >= 0;
psi = gamma1 * (1-indT) + gamma2 * indT;
h = exp(-psi .* abs(deltaT));
if ((isStationary1 == false) | (isStationary2 == false))
    h = h - exp(-(gamma2 * T1 * (isStationary2 == false) ...
        + gamma1 * T2 * (isStationary1 == false)));
end
h = h / (gamma1 + gamma2);
