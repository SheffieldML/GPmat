function [upsilonvp, upsilon] = lfmvpComputeUpsilonMatrix(gamma, sigma2, t1, t2, mode)

% LFMVPCOMPUTEUPSILONMATRIX Upsilon matrix vel. pos. with t1, t2 limits
% FORMAT
% DESC computes a portion of the LFM kernel.
% ARG gamma : Gamma value for system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG mode : operation mode, according to the derivative (mode 0,
% derivative wrt t1, mode 1 derivative wrt t2)
% RETURN upsilon : result of this subcomponent of the kernel for the given values.
%
% SEEALSO : lfmComputeUpsilonMatrix.F, lfmComputeH3.m
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

sigma = sqrt(sigma2);
gridt1 = repmat(t1, 1, length(t2));
gridt2 = repmat(t2', length(t1), 1);
timeGrid = gridt1 - gridt2;

if mode==0
    if nargout > 1
        upsilon = lfmComputeUpsilonMatrix(gamma, sigma2, t1, t2);
        upsilonvp = -gamma*upsilon + (2/(sqrt(pi)*sigma))*exp(-(timeGrid.^2)./sigma2);
    else
        upsilonvp = -gamma*lfmComputeUpsilonMatrix(gamma, sigma2, t1, t2) ...
            + (2/(sqrt(pi)*sigma))*exp(-(timeGrid.^2)./sigma2);
    end
else
    if nargout > 1
        upsilon = lfmComputeUpsilonMatrix(gamma, sigma2, t1, t2);
        upsilonvp = gamma*upsilon - (2/(sqrt(pi)*sigma))*exp(-(timeGrid.^2)./sigma2) ...
            + (2/(sqrt(pi)*sigma))*exp(-gamma*t1)*(exp(-(t2.^2)/sigma2)).';
    else
        upsilonvp = gamma*lfmComputeUpsilonMatrix(gamma, sigma2, t1, t2) ...
            - (2/(sqrt(pi)*sigma))*exp(-(timeGrid.^2)./sigma2) ...
            + (2/(sqrt(pi)*sigma))*exp(-gamma*t1)*(exp(-(t2.^2)/sigma2)).';
    end
end

