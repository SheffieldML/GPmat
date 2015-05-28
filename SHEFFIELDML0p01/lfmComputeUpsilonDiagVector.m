function upsilon = lfmComputeUpsilonDiagVector(gamma, sigma2, t)

% LFMCOMPUTEUPSILONDIAGVECTOR
%
%	Description:
%	


%	Copyright (c) 2010 Mauricio A. Alvarez


upsilon = exp(-gamma*t).*lfmComputeUpsilonVector(-gamma ,sigma2, t);