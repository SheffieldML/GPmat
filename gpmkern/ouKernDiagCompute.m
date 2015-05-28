function k = ouKernDiagCompute(kern, t)

% OUKERNDIAGCOMPUTE Compute diagonal of OU kernel (see ouKernCompute or
% ouKernParamInit for a more detailed description of the OU kernel).
% FORMAT
% DESC computes the diagonal of the kernel matrix for the OU (Ornstein -
% Uhlenbeck kernel given a column vector of inputs. So far the dimension of
% the inputs has to be one.
% ARG kern : the kernel structure for which the kernel matrix is computed.
% ARG t : input data in the form of a column vector.
% RETURN k : a vector of the same size as t containing the diagonal of the
% kernel matrix computed at the given points.
%
% SEEALSO : ouKernParamInit, kernDiagCompute, kernCreate, ouKernCompute
%
% COPYRIGHT : David Luengo, 2009

% KERN


if size(t, 2) > 1
  error('Input can only have one column');
end

c = 0.5*kern.variance/kern.decay;
k = ones(size(t,1), 1);
if (kern.isStationary == false)
    k = k - exp(-2*kern.decay*t);
end
k = c*k;
