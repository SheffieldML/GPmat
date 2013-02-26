function upsilon = lfmComputeUpsilonDiagVector(gamma, sigma2, t)

% LFMCOMPUTEUPSILONDIAGVECTOR 
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% SHEFFIELDML

upsilon = exp(-gamma*t).*lfmComputeUpsilonVector(-gamma ,sigma2, t);