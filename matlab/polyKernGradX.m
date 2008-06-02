function gX = polyKernGradX(kern, X, X2)

% POLYKERNGRADX Gradient of POLY kernel with respect to input locations.
% FORMAT
% DESC computes the gradident of the polynomial
% kernel with respect to the input positions where both the row
% positions and column positions are provided separately.
% ARG kern : kernel structure for which gradients are being
% computed.
% ARG x1 : row locations against which gradients are being computed.
% ARG x2 : column locations against which gradients are being computed.
% RETURN g : the returned gradients. The gradients are returned in
% a matrix which is numData2 x numInputs x numData1. Where numData1 is
% the number of data points in X1, numData2 is the number of data
% points in X2 and numInputs is the number of input
% dimensions in X.
%
% SEEALSO polyKernParamInit, kernGradX, polyKernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN


gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = polyKernGradXpoint(kern, X(i, :), X2);
end
  

function gX = polyKernGradXpoint(kern, x, X2)

% POLYKERNGRADXPOINT Gradient with respect to one point of x.

innerProd = X2*x';  
arg = innerProd*kern.weightVariance + kern.biasVariance;
gX = zeros(size(X2));
for j = 1:size(X2, 2)
  gX(:, j) = kern.degree*kern.weightVariance*kern.variance ...
      *X2(:, j).*arg.^(kern.degree-1);
end
