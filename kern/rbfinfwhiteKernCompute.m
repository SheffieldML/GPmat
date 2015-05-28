function K = rbfinfwhiteKernCompute(kern, t1, t2)

% RBFINFWHITEKERNCOMPUTE Compute the RBF-WHITE kernel (with integration limits
% between minus infinity and infinity) given the parameters, t1 and t2.
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
% SEEALSO : rbfinfwhiteKernParamInit, kernCompute, kernCreate,
% rbfinfwhiteKernDiagCompute, rbfinfwhiteXrbfinfwhiteKernCompute
%
% COPYRIGHT : David Luengo, 2009

% KERN

if nargin < 3
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1-T2;

K = 0.5 * kern.variance * sqrt(kern.inverseWidth/pi) ...
    * exp(- kern.inverseWidth * (deltaT.^2) / 4);
