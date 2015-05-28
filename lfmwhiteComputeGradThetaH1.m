function  gradH = lfmwhiteComputeGradThetaH1(gamma1, gamma2, t1, t2, ...
    gradTheta, isStationary1, isStationary2)

% LFMWHITECOMPUTEGRADTHETAH1 computes a portion of the LFM-WHITE kernel's gradient w.r.t. theta.
% FORMAT
% DESC Helper function for computing part of the gradient of
% the LFM-WHITE kernel w.r.t. to a generic parameter theta (mass, spring or
% damper). Used for obtaining the gradients w.r.t. parameters related to
% the first argument (gamma1).
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG gradTheta : gradient of gamma w.r.t. the generic parameter theta.
% ARG isStationary1: indicates whether the stationary version of the first
% kernel is used (TRUE) or not (FALSE). Set to FALSE by default.
% ARG isStationary2: indicates whether the stationary version of the second
% kernel is used (TRUE) or not (FALSE). Set to FALSE by default.
% RETURN h : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : David Luengo, 2009
%
% SEEALSO : lfmwhiteKernParamInit, lfmwhiteXlfmwhiteKernGradient,
% lfmwhiteComputeH, lfmwhiteComputeGradThetaH2

% KERN


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
