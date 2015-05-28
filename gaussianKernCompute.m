function [K, Kbase, n2]  = gaussianKernCompute(kern, x, x2)

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

if kern.isArd
    sqrtP = sqrt(kern.precisionU);
    sqrtPx = x*sparseDiag(sqrtP);
    if nargin < 3
        n2 = dist2(sqrtPx, sqrtPx);        
    else
        sqrtPx2 = x2*sparseDiag(sqrtP);
        n2 = dist2(sqrtPx, sqrtPx2);        
    end
    Kbase = exp(-0.5*n2);    
else
    if nargin < 3
        n2 = dist2(x, x);        
    else        
        n2 = dist2(x, x2);        
    end
    Kbase = exp(-0.5*kern.precisionU*n2);    
end
K = kern.sigma2Latent*Kbase;    








