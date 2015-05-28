function gT = rbfwhiteKernDiagGradX(kern, t)

% RBFWHITEKERNDIAGGRADX Gradient of RBF-WHITE kernel's diagonal w.r.t. t.
%
%	Description:
%
%	GT = RBFWHITEKERNDIAGGRADX(KERN, T) computes the gradient of the
%	diagonal of the RBF-WHITE kernel matrix with respect to the elements
%	of the input column vector given in t.
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
%	RBFWHITEKERNPARAMINIT, KERNDIAGGRADX, RBFWHITEKERNGRADX


%	Copyright (c) 2009 David Luengo



if size(t, 2) > 1
  error('Input can only have one column');
end

if (kern.isStationary == false)
    gT = 0.5 * kern.variance * kern.inverseWidth / pi ...
        * exp(-kern.inverseWidth * t .* t);
else
    gT = zeros(size(t));
end
