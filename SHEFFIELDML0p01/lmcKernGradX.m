function gX = lmcKernGradX(kern,X, X2)

% LMCKERNGRADX Gradient of LMC kernel with respect to input locations.
%
%	Description:
%
%	G = LMCKERNGRADX(KERN, X, X2) computes the gradient of the LMC
%	kernel with respect to the input positions where both the row
%	positions and column positions are provided separately.
%	 Returns:
%	  G - the returned gradients. The gradients are returned in a matrix
%	   which is numData2 x numInputs x numData1. Where numData1 is the
%	   number of data points in X1, numData2 is the number of data points
%	   in X2 and numInputs is the number of input dimensions in X.
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  X - row locations against which gradients are being computed.
%	  X2 - column locations against which gradients are being computed.
%	
%
%	See also
%	LMCKERNPARAMINIT, KERNGRADX


%	Copyright (c) 2010 Mauricio A. Alvarez


error('Not implemented yet')