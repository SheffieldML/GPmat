function k = ggwhiteKernCompute(kern, x, x2)

% GGWHITEKERNCOMPUTE Compute the GG white kernel given the parameters and X.
% FORMAT
% DESC computes the kernel parameters for
%	the gaussian white gaussian white kernel given inputs associated with rows and
%	columns.
% RETURN k : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG x : the input matrix associated with the rows of the kernel.
% ARG x2 the input matrix associated with the columns of the kernel.
%
% FORMAT
% DESC computes the kernel matrix for the
%	gaussian white kernel given a design matrix of inputs.
%	 Returns:
% RETURN k : the kernel matrix computed at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% ARG X : input data matrix in the form of a design matrix.
%	
% SEEALSO : ggwhiteKernParamInit, kernCompute, kernCreate, ggwhiteKernDiagCompute
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008, 2009

% KERN


if nargin < 3
  x2 = x;
end

Lqr = kern.precisionG;
P = Lqr/2;

if kern.isArd
    sqrtP = sqrt(P);
    sqrtPx = x*sparseDiag(sqrtP);
    sqrtPx2 = x2*sparseDiag(sqrtP);
    n2 = dist2(sqrtPx, sqrtPx2);    
else
    dist = dist2(x, x2);
    n2 = P*dist;
end
factor = kern.sigma2Noise*kern.variance^2;    
Kbase = exp(-0.5*n2);
k = factor*Kbase;

