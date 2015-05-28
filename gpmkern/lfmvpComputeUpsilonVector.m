function [upsilonvp, upsilon] = lfmvpComputeUpsilonVector(gamma, sigma2, t1 ,mode)

% LFMVPCOMPUTEUPSILONVECTOR Upsilon vector for vel. pos. with t1 limit
% FORMAT
% DESC computes a portion of the LFMVP kernel.
% ARG gamma : Gamma value for system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG mode : operation mode, according to the derivative (mode 0,
% derivative wrt t1, mode 1 derivative wrt t2)
% RETURN upsilon : result of this subcomponent of the kernel for the given
% values.
%
% COPYRIGHT : Mauricio A. Alvarez, 2010
%
% SEEALSO : lfmvpComputeUpsilonMatrix.m

% KERN

if nargin<4
    mode = 0;
end

sigma = sqrt(sigma2);
 
if mode==0 
    if nargout > 1
        upsilon = lfmComputeUpsilonVector(gamma, sigma2, t1); 
        upsilonvp = -gamma*upsilon + (2/(sqrt(pi)*sigma))*exp(-(t1.^2)/sigma2);
    else
        upsilonvp = -gamma*lfmComputeUpsilonVector(gamma, sigma2, t1) ...
            + (2/(sqrt(pi)*sigma))*exp(-(t1.^2)/sigma2);
    end
else
    % /~This part is not really needed.
    upsilon = gamma*lfmComputeUpsilonVector(gamma, sigma2, t1) ...
        - (2/(sqrt(pi)*sigma))*exp(-(t1.^2)/sigma2) ...
        + (2/(sqrt(pi)*sigma))*exp(-gamma*t1);
    % ~/
end
