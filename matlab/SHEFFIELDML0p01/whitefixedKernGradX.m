function gX = whitefixedKernGradX(kern, X, X2)

% WHITEFIXEDKERNGRADX Gradient of WHITEFIXED kernel with respect to a point x.
%
%	Description:
%
%	G = WHITEFIXEDKERNGRADX(KERN, X) computes the gradient of the fixed
%	parameter white noise kernel with respect to the input positions.
%	 Returns:
%	  G - the returned gradients. The gradients are returned in a matrix
%	   which is numData x numInputs x numData. Where numData is the
%	   number of data points and numInputs is the number of input
%	   dimensions in X.
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  X - locations against which gradients are being computed.
%
%	G = WHITEFIXEDKERNGRADX(KERN, X1, X2) computes the gradident of the
%	fixed parameter white noise kernel with respect to the input
%	positions where both the row positions and column positions are
%	provided separately.
%	 Returns:
%	  G - the returned gradients. The gradients are returned in a matrix
%	   which is numData2 x numInputs x numData1. Where numData1 is the
%	   number of data points in X1, numData2 is the number of data points
%	   in X2 and numInputs is the number of input dimensions in X.
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  X1 - row locations against which gradients are being computed.
%	  X2 - column locations against which gradients are being computed.
%	
%
%	See also
%	% SEEALSO WHITEFIXEDKERNPARAMINIT, KERNGRADX, WHITEFIXEDKERNDIAGGRADX


%	Copyright (c) 2006 Nathaniel J. King


gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
