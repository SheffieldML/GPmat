function upsilon = lfmComputeUpsilonDiagVector(gamma, sigma2, t)

% LFMCOMPUTEUPSILONDIAGVECTOR 
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% GPMAT

upsilon = exp(-gamma*t).*lfmComputeUpsilonVector(-gamma ,sigma2, t);