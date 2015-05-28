function gX = mlpKernGradX(kern, X, X2)

% MLPKERNGRADX Gradient of MLP kernel with respect to input locations.
% FORMAT
% DESC computes the gradident of the multi-layer perceptron
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
% SEEALSO mlpKernParamInit, kernGradX, mlpKernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009

% KERN


gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = mlpKernGradXpoint(kern, X(i, :), X2);
end
  

function gX = mlpKernGradXpoint(kern, x, X2)

% MLPKERNGRADXPOINT Gradient with respect to one point of x.

innerProd = X2*x';  
numer = innerProd*kern.weightVariance + kern.biasVariance;
vec1 = sum(x.*x, 2)*kern.weightVariance + kern.biasVariance + 1;
vec2 = sum(X2.*X2, 2)*kern.weightVariance + kern.biasVariance + 1;
denom = sqrt(vec2*vec1');
arg = numer./denom;
gX = zeros(size(X2));
twooverpi = 2/pi;
for j = 1:size(X2, 2)
  gX(:, j)=X2(:, j)./denom - vec2.*x(:, j).*numer./denom.^3;
  gX(:, j) = twooverpi*kern.weightVariance*kern.variance*gX(:, j)./sqrt(1-arg.*arg);
end

