function  gradH = lfmwhiteComputeGradThetaH1(gamma1, gamma2, t1, t2, ...
    gradTheta, isStationary1, isStationary2)

% LFMWHITECOMPUTEGRADTHETAH1 computes a portion of the LFM-WHITE kernel's gradient w.r.t. theta.
%
%	Description:
%
%	H = LFMWHITECOMPUTEGRADTHETAH1(GAMMA1, GAMMA2, T1, T2, GRADTHETA,
%	ISSTATIONARY1, ISSTATIONARY2) Helper function for computing part of
%	the gradient of the LFM-WHITE kernel w.r.t. to a generic parameter
%	theta (mass, spring or damper). Used for obtaining the gradients
%	w.r.t. parameters related to the first argument (gamma1).
%	 Returns:
%	  H - result of this subcomponent of the kernel for the given
%	   values.
%	 Arguments:
%	  GAMMA1 - Gamma value for first system.
%	  GAMMA2 - Gamma value for second system.
%	  T1 - first time input (number of time points x 1).
%	  T2 - second time input (number of time points x 1).
%	  GRADTHETA - gradient of gamma w.r.t. the generic parameter theta.
%	  ISSTATIONARY1 - indicates whether the stationary version of the
%	   first kernel is used (TRUE) or not (FALSE). Set to FALSE by
%	   default.
%	  ISSTATIONARY2 - indicates whether the stationary version of the
%	   second kernel is used (TRUE) or not (FALSE). Set to FALSE by
%	   default.
%	
%	lfmwhiteComputeH, lfmwhiteComputeGradThetaH2
%
%	See also
%	LFMWHITEKERNPARAMINIT, LFMWHITEXLFMWHITEKERNGRADIENT, 


%	Copyright (c) 2009 David Luengo



if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

if nargin < 6
    isStationary1 = false;
end
if nargin < 7
    isStationary2 = false;
end

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1 - T2;
indT = deltaT >= 0;
psi = gamma1 * (1-indT) + gamma2 * indT;
gradH = -(lfmwhiteComputeH(gamma1, gamma2, t1, t2, isStationary1, isStationary2) ...
    + abs(deltaT) .* exp(-psi .* abs(deltaT)) .* (1-indT));
if ((isStationary1 == false) | (isStationary2 == false))
    gradH = gradH + T2 .* exp(-(gamma2 * T1 * double(isStationary2 == false) ...
        + gamma1 * T2 * double(isStationary2 == false)));
end
gradH = gradH * gradTheta / (gamma1 + gamma2);
