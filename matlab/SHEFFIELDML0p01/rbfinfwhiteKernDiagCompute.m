function k = rbfinfwhiteKernDiagCompute(kern, t)

% RBFINFWHITEKERNDIAGCOMPUTE Compute diagonal of RBF-WHITE kernel (with
%
%	Description:
%	integration limits between minus infinity and infinity).
%
%	K = RBFINFWHITEKERNDIAGCOMPUTE(KERN, T) computes the diagonal of the
%	kernel matrix for the RBF-WHITE kernel given a column vector of
%	inputs.
%	 Returns:
%	  K - a vector of the same size as t containing the diagonal of the
%	   kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the kernel matrix is
%	   computed.
%	  T - input data in the form of a column vector.
%	rbfinfwhiteKernCompute
%	
%
%	See also
%	RBFINFWHITEKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, 


%	Copyright (c) 2009 David Luengo



if size(t, 2) > 1
  error('Input can only have one column');
end

k = (0.5 * kern.variance * sqrt(kern.inverseWidth/pi)) * ones(size(t));
