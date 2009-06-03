function [K, L]  = gaussianKernCompute(kern, x, x2)

% GAUSSIANKERNCOMPUTE Compute the Gaussian kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for the Gaussian kernel given
% inputs associated with rows and columns.
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : the input matrix associated with the rows of the kernel.
% ARG X2 : the input matrix associated with the columns of the kernel.
%
% FORMAT
% DESC computes the kernel matrix for the Gaussian kernel given a design matrix of inputs.
% RETURN K : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : input data matrix in the form of a design matrix.
%	
% SEEALSO : gaussianKernParamInit, kernCompute, kernCreate, gaussianKernDiagCompute
% 
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS: Mauricio Alvarez, 2009

% KERN

L = sqrt(kern.precision_u);
Lx = x*diag(L);

if nargin < 3  
  n2 = dist2(Lx, Lx);
  K = kern.sigma2_u*exp(-0.5*n2);
else
  Lx2 = x2*diag(L);
  n2 = dist2(Lx, Lx2);
  K = kern.sigma2_u*exp(-0.5*n2);
end

if isfield(kern, 'isNormalised') && ~isempty(kern.isNormalised)
    if kern.isNormalised
        detL = prod(kern.precision_u);
        K = sqrt(detL)*K;        
    end
end








