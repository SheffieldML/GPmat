function gT = ouKernDiagGradX(kern, t)

% OUKERNDIAGGRADX Gradient of OU kernel's diagonal with respect to t (see
%
%	Description:
%	ouKernCompute or ouKernParamInit for a more detailed description of the
%	OU kernel).
%
%	GT = OUKERNDIAGGRADX(KERN, T) computes the gradient of the diagonal
%	of the Ornstein-Uhlenbeck kernel matrix with respect to the elements
%	of the column vector in t.
%	 Returns:
%	  GT - the gradients of the diagonal with respect to each element of
%	   t. The returned vector has the same dimensions as t.
%	 Arguments:
%	  KERN - the kernel structure for which gradients are being
%	   computed.
%	  T - the input data in the form of a design matrix.
%	
%
%	See also
%	OUKERNPARAMINIT, KERNDIAGGRADX, OUKERNGRADX


%	Copyright (c) 2009 David Luengo



if (kern.isStationary == true)
    gT = zeros(size(t));
else
    gT = kern.variance*exp(-2*kern.decay*t);
end
