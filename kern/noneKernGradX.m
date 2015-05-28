function gX = noneKernGradX(kern, x, x2)

% NONEKERNGRADX Gradient of NONE kernel with respect to a point x.
% FORMAT
% DESC computes the gradient of the dummy kernel function
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
% DESC computes the gradident of the dummy kernel function
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
% SEEALSO noneKernParamInit, kernGradX, noneKernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2008

% KERN

if nargin < 3,
    %covGrad = X2;
    x2 = x;
end
gX = zeros(size(x2, 1), size(x2, 2), size(x, 1));
