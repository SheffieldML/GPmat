function dUpsilon = lfmavGradientUpsilonMatrix(gamma, sigma2, t1, ...
    t2, mode, upsilon)

% LFMAVGRADIENTUPSILONMATRIX Gradient upsilon matrix accel. vel.
% FORMAT
% DESC computes the gradient of a portion of the LFMAV kernel.
% ARG gamma : Gamma value for system.
% ARG sigma2 : length scale of latent process.
% ARG t1 : first time input (number of time points x 1).
% ARG t2 : second time input (number of time points x 1).
% ARG upsilon : precomputation of the upsilon matrix.
% ARG mode : operation mode, according to the derivative (mode 0,
% derivative wrt t1, mode 1 derivative wrt t2)
% RETURN upsilon : result of this subcomponent of the kernel for the given values.
%
% COPYRIGHT : Mauricio Alvarez, 2010
%
% SEEALSO : lfmavComputeUpsilonMatrix.m

% KERN

sigma = sqrt(sigma2);

if nargin<6
    upsilon = lfmvpComputeUpsilonMatrix(gamma, sigma2, t1, t2, 1 - mode);
    if nargin <5
        mode =0;
    end
end

dUpsi = lfmvpGradientUpsilonMatrix(gamma, sigma2, t1, t2, 1 - mode);
gridt1 = repmat(t1, 1, length(t2));
gridt2 = repmat(t2', length(t1), 1);
timeGrid = gridt1 - gridt2;

if mode ==0
    dUpsilon = 2*gamma*upsilon + (gamma^2)*dUpsi - (4/(sqrt(pi)*sigma^3))* ...
        timeGrid.*exp(-(timeGrid.^2)./sigma2);
else
    dUpsilon = 2*gamma*upsilon + (gamma^2)*dUpsi + (4/(sqrt(pi)*sigma^3))* ...
        timeGrid.*exp(-(timeGrid.^2)./sigma2) - (2/(sqrt(pi)*sigma))*...
        (exp(-gamma*t1).*(1-gamma*t1))*((gamma - 2*t2/sigma2).*exp(-(t2.^2)/sigma2)).' ...
        - (2*gamma/(sqrt(pi)*sigma))*exp(-gamma*t1)*(exp(-(t2.^2)/sigma2)).';        
end
