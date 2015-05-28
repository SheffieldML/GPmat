function [upsilonaa, upsilonap] = lfmaaComputeUpsilonMatrix(gamma, sigma2, t1, t2, mode)

% LFMAACOMPUTEUPSILONMATRIX Upsilon matrix acce. accel. with t1, t2 limits
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
% COPYRIGHT : Mauricio Alvarez, 2010

% KERN

sigma = sqrt(sigma2);
gridt1 = repmat(t1, 1, length(t2));
gridt2 = repmat(t2', length(t1), 1);
timeGrid = gridt1 - gridt2;

if mode==0
    if nargout > 1
        upsilonap = lfmapComputeUpsilonMatrix(gamma, sigma2, t1, t2,1);
        upsilonaa = gamma^2*upsilonap + (2/(sqrt(pi)*sigma))*exp(-(timeGrid.^2)./sigma2).* ...
            ((gamma + 2*timeGrid/sigma2).*(2/sigma2 - (2*timeGrid/sigma2).^2) + 8*timeGrid/sigma^4);
    else
        upsilonaa = gamma^2*lfmapComputeUpsilonMatrix(gamma, sigma2, t1, t2,1) ...
            + (2/(sqrt(pi)*sigma))*exp(-(timeGrid.^2)./sigma2).* ...
            ((gamma + 2*timeGrid/sigma2).*(2/sigma2 - (2*timeGrid/sigma2).^2) + 8*timeGrid/sigma^4);
    end
else
    if nargout > 1
        upsilonap = lfmapComputeUpsilonMatrix(gamma, sigma2, t1, t2, 0);
        upsilonaa = gamma^2*upsilonap + (2/(sqrt(pi)*sigma))*exp(-(timeGrid.^2)./sigma2).* ...
            ((gamma + 2*timeGrid/sigma2).*(2/sigma2 - (2*timeGrid/sigma2).^2) + 8*timeGrid/sigma^4) ...
            + (2*gamma^2/(sqrt(pi)*sigma))*exp(-gamma*t1)*((gamma - 2*t2/sigma2).*exp(-(t2.^2)/sigma2)).';

    else
        upsilonaa = gamma^2*lfmapComputeUpsilonMatrix(gamma, sigma2, t1, t2, 0) ...
            + (2/(sqrt(pi)*sigma))*exp(-(timeGrid.^2)./sigma2).* ...
            ((gamma + 2*timeGrid/sigma2).*(2/sigma2 - (2*timeGrid/sigma2).^2) + 8*timeGrid/sigma^4) ...
            + (2*gamma^2/(sqrt(pi)*sigma))*exp(-gamma*t1)*((gamma - 2*t2/sigma2).*exp(-(t2.^2)/sigma2)).';
    end
end
