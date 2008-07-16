function g = ivmCovarianceGradient(invK, m)

% IVMCOVARIANCEGRADIENT The gradient of the likelihood approximation wrt the covariance.
% FORMAT
% DESC returns the gradient of the log likelihood approximation
% with respect to the covariance function (i.e. kernel) evaluated
% at the active points.
% ARG invK : the inverse of the covariance function evaluated at
% the active points.
% ARG m : the difference between the target and the mean value.
% RETURN g : the gradient of the log likelihood approximation with
% respect to the covariance function elements.
%
% SEEALSO : ivmCreate
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% IVM

invKm = invK*m;

g = -invK + invKm*invKm';
g= g*.5;
