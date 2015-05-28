function k = kernCorrelation(kern, x, x2)

% KERNCORRELATION Compute the correlation matrix kernel given the parameters and X.
% FORMAT
% DESC computes a correlation matrix for the given kernel type given an
% input data matrix.
% ARG kern : kernel structure to be computed.
% ARG X : input data matrix (rows are data points) to the kernel computation.
% RETURN K : computed elements of the correlation matrix.
%
% FORMAT 
% DESC computes a correlation matrix for the given kernel type given 
% two input data matrices, one for the rows and one for the columns.
% ARG kern : kernel structure to be computed.
% ARG X : first input matrix to the kernel computation (forms the rows of the kernel).
% ARG X2 : second input matrix to the kernel computation (forms the columns of the kernel).
% RETURN K : computed elements of the correlation matrix.
%
% SEEALSO : kernCreate, kernDiagCompute, kernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2007
  
% KERN

fhandle = str2func([kern.type 'KernCompute']);
if nargin < 3
  k = fhandle(kern, x);
  if ~kern.isStationary
    d = 1./sqrt(kernDiagCompute(kern, x));
    k = sparseDiag(d)*k*sparseDiag(d);
  end
else
  k = fhandle(kern, x, x2);
  if ~kern.isStationary
    k = sparseDiag(1./sqrt(kernDiagCompute(kern, x)))...
        *k*sparseDiag(1./sqrt(kernDiagCompute(kern, x2)));
  end
end
