function gT = rbfinfwhiteKernDiagGradX(kern, t)

% RBFINFWHITEKERNDIAGGRADX Gradient of RBF-WHITE kernel's (with integration
%
%	Description:
%	limits between minus infinity and infinity) diagonal w.r.t. t.
%
%	GT = RBFINFWHITEKERNDIAGGRADX(KERN, T) computes the gradient of the
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
%	RBFINFWHITEKERNPARAMINIT, KERNDIAGGRADX, RBFINFWHITEKERNGRADX


%	Copyright (c) 2009 David Luengo



if size(t, 2) > 1
  error('Input can only have one column');
end

gT = zeros(size(t));
