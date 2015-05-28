function gX = whiteblockKernGradX(kern, X, X2)

% WHITEBLOCKKERNGRADX Gradient of WHITEBLOCK kernel wrt input locations.
%
%	Description:
%
%	G = WHITEBLOCKKERNGRADX(KERN, X1, X2) computes the gradident of the
%	white noise block kernel with respect to the input positions where
%	both the row positions and column positions are provided separately.
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
%	% SEEALSO WHITEBLOCKKERNPARAMINIT, KERNGRADX, WHITEBLOCKKERNDIAGGRADX


%	Copyright (c) 2010 Mauricio A. Alvarez


if nargin<3
    X2 = X;
end

gX = cell(1, kern.nout);
for i=1:kern.nout
    gX{i} = zeros(size(X2, 1), size(X2, 2), size(X, 1));
end