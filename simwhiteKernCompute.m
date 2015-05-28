function K = simwhiteKernCompute(kern, t1, t2)

% SIMWHITEKERNCOMPUTE Compute the SIM-WHITE kernel given the parameters, t1
% and t2. This kernel is the result of using a white noise latenf function
% to force the SIM kernel.
% FORMAT
% DESC computes the kernel matrix for the SIM-White (Single Input Motif
% - White) kernel given inputs associated with rows and columns.
% ARG kern : the kernel structure for which the kernel matrix is computed.
% ARG t1 : the input column vector associated with the rows of the kernel.
% ARG t2 : the input column vector associated with the columns of the kernel.
% RETURN K : the kernel matrix computed at the given points.
%
% FORMAT
% DESC computes the kernel matrix for the SIM-White (Single Input Motif - White)
% kernel given a design matrix of inputs.
% ARG kern : the kernel structure for which the kernel matrix is computed.
% ARG t1 : input data in the form of a column vector.
% RETURN K : the kernel matrix computed at the given points.
%
% SEEALSO : simwhiteKernParamInit, kernCompute, kernCreate, simwhiteKernDiagCompute
%
% COPYRIGHT : David Luengo, 2009

% KERN

if nargin < 3
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

c = kern.variance*(kern.sensitivity^2)/(2*kern.decay);
T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
K = exp(-kern.decay*abs(T1-T2));
if (kern.isStationary == false)
    K = K - exp(-kern.decay*(T1+T2));
end
K = c*K;
