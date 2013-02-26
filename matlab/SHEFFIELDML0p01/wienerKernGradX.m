function gX = wienerKernGradX(kern, x, x2)

% WIENERKERNGRADX Gradient of WIENER kernel with respect to a point x.
%
%	Description:
%
%	G = WIENERKERNGRADX(KERN, X) computes the gradient of the wiener
%	kernel with respect to the input positions.
%	 Returns:
%	  G - the returned gradients. The gradients are returned in a matrix
%	   which is numData x numInputs x numData. Where numData is the
%	   number of data points and numInputs is the number of input
%	   dimensions in X.
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  X - locations against which gradients are being computed.
%
%	G = WIENERKERNGRADX(KERN, X1, X2) computes the gradident of the
%	wiener kernel with respect to the input positions where both the row
%	positions and column positions are provided separately.
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
%	% SEEALSO WIENERKERNPARAMINIT, KERNGRADX, WIENERKERNDIAGGRADX


%	Copyright (c) 2009 Neil D. Lawrence


error(['Unable to compute gradient with respect to input variables for Wiener kernel.']);