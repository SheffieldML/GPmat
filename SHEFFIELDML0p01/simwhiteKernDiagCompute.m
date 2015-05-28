function K = simwhiteKernDiagCompute(kern, t)

% SIMWHITEKERNDIAGCOMPUTE Compute the diagonal of the SIM-WHITE kernel.
%
%	Description:
%
%	K = SIMWHITEKERNDIAGCOMPUTE(KERN, T) computes the diagonal of the
%	kernel matrix for the SIM-White (Single Input Motif - White) kernel
%	given a column vector of inputs.
%	 Returns:
%	  K - a vector of the same size as t containing the diagonal of the
%	   kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the kernel matrix is
%	   computed.
%	  T - input data in the form of a column vector.
%	
%
%	See also
%	SIMWHITEKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, SIMWHITEKERNCOMPUTE


%	Copyright (c) 2009 David Luengo



if size(t, 2) > 1
  error('Input can only have one column');
end

c = kern.variance * (kern.sensitivity^2) / (2*kern.decay);
K = ones(size(t));
if (kern.isStationary == false)
    K = K - exp(-2*kern.decay*t);
end
K = c*K;
