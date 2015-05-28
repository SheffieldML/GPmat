function gX = polyardKernGradX(kern, X, X2)

% POLYARDKERNGRADX Gradient of POLYARD kernel with respect to input locations.
%
%	Description:
%
%	G = POLYARDKERNGRADX(KERN, X1, X2) computes the gradident of the
%	automatic relevance determination polynomial kernel with respect to
%	the input positions where both the row positions and column
%	positions are provided separately.
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
%	% SEEALSO POLYARDKERNPARAMINIT, KERNGRADX, POLYARDKERNDIAGGRADX


%	Copyright (c) 2005, 2006 Neil D. Lawrence



gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = polyardKernGradXpoint(kern, X(i, :), X2);
end
  
function gX = polyardKernGradXpoint(kern, x, X2)

% POLYARDKERNGRADXPOINT Gradient with respect to one point of x.

scales = sparse(diag(kern.inputScales));
xScaled = x*scales;
X2Scaled = X2*scales;
innerProd = X2Scaled*x';
arg = innerProd*kern.weightVariance + kern.biasVariance;
gX = zeros(size(X2));
for j = 1:size(X2, 2)
  gX(:, j) = kern.degree*kern.weightVariance*kern.inputScales(j)*kern.variance*X2(:, j).*arg.^(kern.degree-1);
end
