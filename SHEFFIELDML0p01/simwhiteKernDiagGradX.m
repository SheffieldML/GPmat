function gT = simwhiteKernDiagGradX(kern, t)

% SIMWHITEKERNDIAGGRADX Gradient of SIM-WHITE kernel's diagonal w.r.t. t.
%
%	Description:
%
%	GT = SIMWHITEKERNDIAGGRADX(KERN, T) computes the gradient of the
%	diagonal of the SIM-White (Single Input Motif - White) kernel matrix
%	with respect to the elements of the input column vector given in t.
%	 Returns:
%	  GT - the gradients of the diagonal with respect to each element of
%	   t. The returned matrix has the same dimensions as t.
%	 Arguments:
%	  KERN - the kernel structure for which gradients are being
%	   computed.
%	  T - the input data in the form of a design matrix.
%	
%
%	See also
%	SIMWHITEKERNPARAMINIT, KERNDIAGGRADX, SIMWHITEKERNGRADX


%	Copyright (c) 2009 David Luengo



if (kern.isStationary == true)
    gT = zeros(size(t));
else
    gT = kern.variance * (kern.sensitivity^2) * exp(-2*kern.decay*t);
end
