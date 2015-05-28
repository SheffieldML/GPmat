function gX = translateKernGradX(kern, varargin)

% TRANSLATEKERNGRADX Gradient of TRANSLATE kernel with respect to a point x.
%
%	Description:
%
%	G = TRANSLATEKERNGRADX(KERN, X) computes the gradient of the input
%	space translation kernel with respect to the input positions.
%	 Returns:
%	  G - the returned gradients. The gradients are returned in a matrix
%	   which is numData x numInputs x numData. Where numData is the
%	   number of data points and numInputs is the number of input
%	   dimensions in X.
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  X - locations against which gradients are being computed.
%
%	G = TRANSLATEKERNGRADX(KERN, X1, X2) computes the gradident of the
%	input space translation kernel with respect to the input positions
%	where both the row positions and column positions are provided
%	separately.
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
%	% SEEALSO TRANSLATEKERNPARAMINIT, KERNGRADX, CMPNDKERNGRADX, TRANSLATEKERNDIAGGRADX


%	Copyright (c) 2007 Neil D. Lawrence


for i = 1:length(varargin)
  varargin{i} = varargin{i} - repmat(kern.centre, size(varargin{i}, 1), 1);
end
gX = cmpndKernGradX(kern, varargin{:});
