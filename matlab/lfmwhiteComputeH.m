function  h = lfmwhiteComputeH(gamma1, gamma2, t1, t2, isStationary1, isStationary2)

% LFMWHITECOMPUTEH Helper function for computing part of the LFM-WHITE
% kernel.
% FORMAT
% DESC computes a portion of the LFM-WHITE kernel.
% ARG gamma1 : Gamma value for first system.
% ARG gamma2 : Gamma value for second system.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG isStationary1: indicates whether the stationary version of the first
% kernel is used (1) or not (0).
% ARG isStationary2: indicates whether the stationary version of the second
% kernel is used (1) or not (0).
% RETURN h : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : David Luengo, 2009
%
% SEEALSO : lfmwhiteKernParamInit, lfmwhiteXlfmwhiteKernCompute

% KERN


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
