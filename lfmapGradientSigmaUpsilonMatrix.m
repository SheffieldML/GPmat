function dUpsilonS = lfmapGradientSigmaUpsilonMatrix(gamma, sigma2, ...
    t1, t2, mode)

% LFMAPGRADIENTSIGMAUPSILONMATRIX Gradient of upsilon matrix ap wrt sigma
% FORMAT
% DESC computes the gradient of a portion of the LFMAP kernel.
% ARG gamma : Gamma value for system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG mode : operation mode, according to the derivative (mode 0,
% derivative wrt t1, mode 1 derivative wrt t2)
% RETURN upsilon : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio Alvarez, 2010
%
% SEEALSO : lfmvpComputeUpsilonMatrix.m

% KERN

gridt1 = repmat(t1, 1, length(t2));
gridt2 = repmat(t2', length(t1), 1);
timeGrid = gridt1 - gridt2;

dUpsilon = lfmGradientSigmaUpsilonMatrix(gamma, sigma2, t1, t2);

if mode == 0
    dUpsilonS = (gamma^2)*dUpsilon + (2/(sqrt(pi)*sigma2))* ...
        exp(-(timeGrid.^2)./sigma2).*(gamma + 2*timeGrid/sigma2).*(1-2*(timeGrid.^2)/sigma2) ...
        + (8/(sqrt(pi)*sigma2^2))*timeGrid.*exp(-(timeGrid.^2)./sigma2);
else
    dUpsilonS = (gamma^2)*dUpsilon + (2/(sqrt(pi)*sigma2))* ...
        exp(-(timeGrid.^2)./sigma2).*(gamma + 2*timeGrid/sigma2).*(1-2*(timeGrid.^2)/sigma2) ...
        + (8/(sqrt(pi)*sigma2^2))*timeGrid.*exp(-(timeGrid.^2)./sigma2) ...
        - (2/(sqrt(pi)*sigma2))*exp(-gamma*t1)*((gamma-2*t2/sigma2).*(1-2*t2.^2/sigma2).*exp(-t2.^2/sigma2)).' ...
        + (8/(sqrt(pi)*sigma2^2))*exp(-gamma*t1)*(t2.*exp(-t2.^2/sigma2)).';
end
