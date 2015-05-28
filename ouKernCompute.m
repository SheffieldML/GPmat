function [K, sk] = ouKernCompute(kern, t1, t2)

% OUKERNCOMPUTE Compute the Ornstein-Uhlenbeck (OU) kernel arising from the
% stochastic differential equation
%
%    df(t) = -theta*(f(t)-mu)dt + sigma*dw(t),
%
% where w(t) is the Wiener process, given the parameters (theta and sigma),
% and t.
% FORMAT
% DESC computes the kernel parameters for the Ornstein-Uhlenbeck kernel
% given inputs associated with rows and columns. So far the dimension of
% the inputs has to be one.
% ARG kern : the kernel structure for which the kernel matrix is computed.
% ARG t1 : the input column vector associated with the rows of the kernel.
% ARG t2 : the input column vector associated with columns of the kernel.
% RETURN K : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the Ornstein-Uhlenbeck kernel
% given a column vector of inputs.
% ARG kern : the kernel structure for which the kernel matrix is computed.
% ARG t1 : input data in the form of a column vector.
% RETURN K : the kernel matrix computed at the given points.
%
% SEEALSO : ouKernParamInit, kernCompute, kernCreate, ouKernDiagCompute
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
c = 0.5/kern.decay;
K = exp(-kern.decay*abs(T1-T2));
if (kern.isStationary == false)
    K = K - exp(-kern.decay*(T1+T2));
end
sk = c*K;
K = sk*kern.variance;
