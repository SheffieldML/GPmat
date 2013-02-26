function gT = lfmwhiteKernDiagGradX(kern, t)

% LFMWHITEKERNDIAGGRADX Gradient of LFM-WHITE kernel's diagonal w.r.t. t.
%
%	Description:
%
%	GT = LFMWHITEKERNDIAGGRADX(KERN, T) computes the gradient of the
%	diagonal of the SIM-White (Latent Force Model - White) kernel matrix
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
%	LFMWHITEKERNPARAMINIT, KERNDIAGGRADX, LFMWHITEKERNGRADX


%	Copyright (c) 2009 David Luengo



if size(t, 2) > 1
  error('Input can only have one column');
end

if (kern.isStationary == false)
    c = kern.variance*(kern.sensitivity^2) / (4*(kern.mass^2)*(kern.omega^2));
    gT = c * (2 - exp(-j*2*kern.omega*t) - exp(j*2*kern.omega*t)) ...
        .* exp(-2*kern.alpha*t);
else
    gT = zeros(size(t));
end
