function gX = velotransKernGradX(kern, X, X2)

% VELOTRANSKERNGRADX Gradient of VELOTRANS kernel with respect to a point x.
%
%	Description:
%
%	G = VELOTRANSKERNGRADX(KERN, X) computes the gradient of the
%	velocity translate kernel with respect to the input positions.
%	 Returns:
%	  G - the returned gradients. The gradients are returned in a matrix
%	   which is numData x numInputs x numData. Where numData is the
%	   number of data points and numInputs is the number of input
%	   dimensions in X.
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  X - locations against which gradients are being computed.
%
%	G = VELOTRANSKERNGRADX(KERN, X1, X2) computes the gradident of the
%	velocity translate kernel with respect to the input positions where
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
%	% SEEALSO VELOTRANSKERNPARAMINIT, KERNGRADX, VELOTRANSKERNDIAGGRADX, TRANSLATEKERNGRADX


%	Copyright (c) 2011 Neil D. Lawrence


  t = X(:, end);
  xPass = X(:, 1:end-1) - t*kern.velocity;
  t2 = X2(:, end);
  x2Pass = X2(:, 1:end-1) - t2*kern.velocity;
  gXRecover = cmpndKernGradX(kern, xPass, x2Pass);
  gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
  for i = 1:size(X, 1)
    gX(:, 1:end-1, i) = gXRecover(:, :, i);
    gX(:, end, i) = -gXRecover(:, :, i)*kern.velocity'; % Another gradient
                                                        % guess
  end
end