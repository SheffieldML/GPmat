function gX = rbfardjitKernGradX(kern, X, X2)

% RBFARDJITKERNGRADX Gradient of RBFARDJIT kernel with respect to input locations.
%
%	Description:
%
%	G = RBFARDJITKERNGRADX(KERN, X1, X2) computes the gradident of the
%	automatic relevance determination radial basis function kernel with
%	respect to the input positions where both the row positions and
%	column positions are provided separately.
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
%
%	See also
%	% SEEALSO RBFARD2KERNPARAMINIT, KERNGRADX, RBFARD2KERNDIAGGRADX


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2009 Michalis K. Titsias



gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = rbfardjitKernGradXpoint(kern, X(i, :), X2);
end
  

function gX = rbfardjitKernGradXpoint(kern, x, X2)

% RBFARDJITKERNGRADXPOINT Gradient with respect to one point of x.

scales = sparse(sqrt(diag(kern.inputScales)));
gX = zeros(size(X2));
n2 = dist2(X2*scales, x*scales);
rbfPart = kern.variance*exp(-n2*0.5);
for i = 1:size(x, 2)
  gX(:, i) = kern.inputScales(i)*(X2(:, i) - x(i)).*rbfPart;
end
