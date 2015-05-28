function K = rbfwhiteKernCompute(kern, t1, t2)

% RBFWHITEKERNCOMPUTE Compute the RBF-WHITE kernel given the parameters, t1
% and t2.
% FORMAT
% DESC computes the kernel parameters for the RBF-WHITE kernel given inputs
% associated with rows and columns.
% ARG kern : the kernel structure for which the kernel matrix is computed.
% ARG t1 : the input column vector associated with the rows of the kernel.
% ARG t2 : the input column vector associated with the columns of the kernel.
% RETURN K : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the RBF-WHITE kernel given a design
% matrix of inputs.
% ARG kern : the kernel structure for which the kernel matrix is computed.
% ARG t1 : input data in the form of a column vector.
% RETURN K : the kernel matrix computed at the given points.
%
% SEEALSO : rbfwhiteKernParamInit, kernCompute, kernCreate,
% rbfwhiteKernDiagCompute, rbfwhiteXrbfwhiteKernCompute
%
% COPYRIGHT : David Luengo, 2009

% KERN

if nargin < 3
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

K = rbfwhiteXrbfwhiteKernCompute(kern, kern, t1, t2);

if nargin < 3;
  K = 0.5*(K + K');
end

K = real(K);
