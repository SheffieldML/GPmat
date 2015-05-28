function gX = whiteblockKernGradX(kern, X, X2)

% WHITEBLOCKKERNGRADX Gradient of WHITEBLOCK kernel wrt input locations.
% FORMAT
% DESC computes the gradident of the white noise block kernel with respect 
% to the input positions where both the row positions and column positions 
% are provided separately.
% ARG kern : kernel structure for which gradients are being computed.
% ARG x1 : row locations against which gradients are being computed.
% ARG x2 : column locations against which gradients are being computed.
% RETURN g : the returned gradients. The gradients are returned in
% a matrix which is numData2 x numInputs x numData1. Where numData1 is
% the number of data points in X1, numData2 is the number of data
% points in X2 and numInputs is the number of input
% dimensions in X.
%
% SEEALSO whiteblockKernParamInit, kernGradX, whiteblockKernDiagGradX
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if nargin<3
    X2 = X;
end

gX = cell(1, kern.nout);
for i=1:kern.nout
    gX{i} = zeros(size(X2, 1), size(X2, 2), size(X, 1));
end
