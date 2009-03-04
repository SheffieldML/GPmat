function [K, Pinv]  = gaussianwhiteKernCompute(kern, x, x2)

% GAUSSIANWHITEKERNCOMPUTE Compute the covariance of the output samples 
% when the input is a white noise process and the smoothing kernel is a 
% Gaussian kernel.
% FORMAT
% DESC computes the kernel parameters for the Gaussian white kernel given
% inputs associated with rows and columns.
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : the input matrix associated with the rows of the kernel.
% ARG X2 : the input matrix associated with the columns of the kernel.
%
% FORMAT
% DESC computes the kernel matrix for the Gaussian white kernel given a design matrix of inputs.
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
%	
% SEEALSO : gaussianwhiteKernParamInit, kernCompute, kernCreate, gaussianwhiteKernDiagCompute
% 
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008

% KERN

sqrtP = sqrt(kern.precisionT/2);
Pinv = 2./kern.precisionT;
Px = x*sparseDiag(sqrtP);

if nargin < 3  
  n2 = dist2(Px, Px);
  factor = kern.sigma2Noise/((2*pi)^(kern.inputDimension/2)*sqrt(prod(Pinv))); 
  K = factor*exp(-0.5*n2);
else
  Px2 = x2*sparseDiag(sqrtP);
  n2 = dist2(Px, Px2);
  factor = kern.sigma2Noise/((2*pi)^(kern.inputDimension/2)*sqrt(prod(Pinv))); 
  K = factor*exp(-0.5*n2);
end










