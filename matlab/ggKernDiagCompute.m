function K = ggKernDiagCompute(kern, x)

% GGKERNDIAGCOMPUTE Compute diagonal of GG kernel.
% FORMAT
% DESC computes the diagonal of the kernel
%	matrix for the gaussian gaussian kernel given a design matrix of
%	inputs.
% RETURN  K : a vector containing the diagonal of the kernel matrix computed
%	   at the given points.
% ARG kern : the kernel structure for which the matrix is computed.
% X : input data matrix in the form of a design matrix.
%	
% SEEALSO : ggKernParamInit, kernDiagCompute, kernCreate, ggKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008

% KERN



Ank = kern.precision_y;
Bk = kern.precision_u;
Ankinv = 1./Ank;
Bkinv = 1./Bk;
detBkinv = prod(Bkinv);
P = 2*Ankinv + Bkinv;
ldet = prod(P);
K = repmat(kern.sigma2_y^2*kern.sigma2_u*sqrt((detBkinv)/ldet), ...
    size(x,1),1);
