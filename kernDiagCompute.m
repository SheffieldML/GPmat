function k = kernDiagCompute(kern, x)

% KERNDIAGCOMPUTE Compute the kernel given the parameters and X.
% FORMAT
% DESC computes the diagonal of a kernel matrix for the given kernel.
% ARG kern : kernel structure to be computed.
% ARG X : input data matrix (rows are data points) to the kernel computation.
% RETURN K : vector containing computed diagonal elements of the
% kernel structure.
%
% SEEALSO : kernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004, 2005, 2006

% KERN

fhandle = str2func([kern.type 'KernDiagCompute']);
k = fhandle(kern, x);
