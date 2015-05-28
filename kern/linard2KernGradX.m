function gX = linard2KernGradX(kern, X, X2)

% LINARD2KERNGRADX Gradient of LINARD2 kernel with respect to input locations.
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
% SEEALSO linard2KernParamInit, kernGradX, linard2KernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% COPYRIGHT : Michalis K. Titsias, 2009
% MODIFICATIONS: Andreas C. Damianou, 2012

% KERN


scales = sparse(diag(kern.inputScales));

%%%-- Allow for 1-D inputs (Q=1)
if isscalar(scales)
    scales = full(scales);
end
%%%--

X2 = X2*scales;

gX = repmat(X2, [1 1 size(X, 1)]);
