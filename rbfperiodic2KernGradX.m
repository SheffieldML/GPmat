function gX = rbfperiodic2KernGradX(kern, X, X2)


% RBFPERIODIC2KERNGRADX Gradient of RBFPERIODIC2 kernel with respect to a point x.
% FORMAT
% DESC computes the gradient of the RBF periodic covariance with variying period
% kernel with respect to the input positions. 
% ARG kern : kernel structure for which gradients are being
% computed.
% ARG x : locations against which gradients are being computed.
% RETURN g : the returned gradients. The gradients are returned in
% a matrix which is numData x numInputs x numData. Where numData is
% the number of data points and numInputs is the number of input
% dimensions in X.
%
% FORMAT
% DESC computes the gradident of the RBF periodic covariance with variying period
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
% SEEALSO rbfperiodic2KernParamInit, kernGradX, rbfperiodic2KernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2007, 2009
%
% MODIFICATIONS : Andreas C. Damianou, 2011
%
% MODIFICATIONS : Michalis K. Titsias, 2011

% KERN



gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = rbfperiodic2KernGradXpoint(kern, X(i, :), X2);
end
  

function gX = rbfperiodic2KernGradXpoint(kern, x, X2)

% PERIODICKERNGRADXPOINT Gradient with respect to one point of x.

factor = kern.factor;
gX = zeros(size(X2));
arg = factor*(X2 - x)/2;
sinarg = sin(arg);
n2 = sinarg.*sinarg;
wi2 = (2 .* kern.inverseWidth);
rbfPart = kern.variance*exp(-n2*wi2);

gX = factor*kern.inverseWidth*cos(arg).*sinarg.*2.*rbfPart;
