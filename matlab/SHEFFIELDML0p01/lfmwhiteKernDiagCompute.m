function k = lfmwhiteKernDiagCompute(kern, t)

% LFMWHITEKERNDIAGCOMPUTE Compute diagonal of LFM-WHITE kernel.
%
%	Description:
%
%	K = LFMWHITEKERNDIAGCOMPUTE(KERN, T) computes the diagonal of the
%	kernel matrix for the LFM-White (Latent Force Model - White) kernel
%	given a column vector of inputs.
%	 Returns:
%	  K - a vector of the same size as t containing the diagonal of the
%	   kernel matrix computed at the given points.
%	 Arguments:
%	  KERN - the kernel structure for which the kernel matrix is
%	   computed.
%	  T - input data in the form of a column vector.
%	lfmwhiteKernCompute
%	
%
%	See also
%	LFMWHITEKERNPARAMINIT, KERNDIAGCOMPUTE, KERNCREATE, 


%	Copyright (c) 2009 David Luengo



if size(t, 2) > 1
  error('Input can only have one column');
end

gamma2 = kern.alpha - j*kern.omega;
c = kern.variance * (kern.sensitivity^2) / (4*(kern.mass^2)*(kern.omega^2));
k = 1/kern.alpha - kern.alpha/((kern.alpha^2)+(kern.omega^2));
if (kern.isStationary == false)
    k = k - exp(-2*kern.alpha*t) / kern.alpha ...
        + (kern.gamma * exp(-2*gamma2*t) + gamma2 * exp(-2*kern.gamma*t)) ...
        / (2*((kern.alpha^2)+(kern.omega^2)));
end
k = c*real(k);
