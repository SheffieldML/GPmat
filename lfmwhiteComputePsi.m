function psi = lfmwhiteComputePsi(gamma, invWidth, t1, t2, isStationary)

% LFMWHITECOMPUTEPSI Helper function for comptuing part of the LFM-WHITE
% kernel.
% FORMAT
% DESC computes a portion of the LFM-WHITE kernel.
% ARG gamma : Gamma value for the system.
% ARG invWidth : inverse width of the RBF smoothing process.
% ARG t1 : first time input (number of time points of x1).
% ARG t2 : second time input (number of time points of x2).
% ARG isStationary : indicates whether the stationary or non-stationary
% version of the function is used.
% RETURN psi : result of this subcomponent of the kernel for the given
% values.
%
% COPYRIGHT : David Luengo, 2009
%
% SEEALSO : lfmwhiteKernParamInit, lfmwhiteXlfmwhiteKernCompute, W

% KERN


if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

% Time matrices and difference
T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1-T2;
indT = double(T1 >= T2);

% Auxiliary functions and computation of psi (stationary case)

varphiT1T2 = gamma * deltaT .* indT + 0.5 * invWidth * (deltaT.^2) .* (1-indT);
zT1T2 = sqrt(0.5*invWidth) * (-deltaT .* (1-indT) + gamma/invWidth);
psi = exp(-varphiT1T2) .* wofzHui(j*zT1T2);

% Additional auxiliary functions and terms of psi (non-stationary case)

if (isStationary == false)
    varphiT10 = gamma * T1;
    varphi0T2 = 0.5 * invWidth * (T2.^2);
    z0T2 = sqrt(0.5*invWidth) * (T2 + gamma/invWidth);
    psi = psi - exp(-(varphiT10 + varphi0T2)) .* wofzhui(j*z0T2);
end
