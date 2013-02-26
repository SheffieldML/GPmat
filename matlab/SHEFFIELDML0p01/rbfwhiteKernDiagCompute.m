function k = rbfwhiteKernDiagCompute(kern, t)

% RBFWHITEKERNDIAGCOMPUTE Compute diagonal of RBF-WHITE kernel.
%
%	Description:
%
%	K = RBFWHITEKERNDIAGCOMPUTE(KERN, T) computes the diagonal of the
%	kernel matrix for the RBF-WHITE kernel given a column vector of
%	inputs.
%	 Returns:
%	  K - a vector of the same size as t containing the diagonal of the
%	   kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the kernel matrix is
%	   computed.
%	  T - input data in the form of a column vector.
%	rbfwhiteKernCompute
%	
%
%	See also
%	RBFWHITEKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, 


%	Copyright (c) 2009 David Luengo



if size(t, 2) > 1
  error('Input can only have one column');
end

if (kern.isStationary == false)
    k = kern.variance/sqrt(8*pi) * erf(kern.inverseWidth*t/sqrt(2));
else
    k = zeros(size(t));
end
