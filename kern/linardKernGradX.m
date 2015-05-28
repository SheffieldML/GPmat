function gX = linardKernGradX(kern, X, X2)

% LINARDKERNGRADX Gradient of LINARD kernel with respect to input locations.
% FORMAT
% DESC computes the gradident of the automatic relevance determination linear
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
% SEEALSO linardKernParamInit, kernGradX, linardKernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


scales = sparse(diag(kern.inputScales));
X2 = X2*scales;

gX = repmat(kern.variance.*X2, [1 1 size(X, 1)]);
