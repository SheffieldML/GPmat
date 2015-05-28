function upsilon = lfmComputeUpsilonDiagVector(gamma, sigma2, t)

% LFMCOMPUTEUPSILONDIAGVECTOR 
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

upsilon = exp(-gamma*t).*lfmComputeUpsilonVector(-gamma ,sigma2, t);
