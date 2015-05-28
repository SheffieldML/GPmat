function gX = translateKernGradX(kern, varargin)

% TRANSLATEKERNGRADX Gradient of TRANSLATE kernel with respect to a point x.
% FORMAT
% DESC computes the gradient of the input space translation
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
% DESC computes the gradident of the input space translation
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
% SEEALSO translateKernParamInit, kernGradX, cmpndKernGradX, translateKernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN

for i = 1:length(varargin)
  varargin{i} = varargin{i} - repmat(kern.centre, size(varargin{i}, 1), 1);
end
gX = cmpndKernGradX(kern, varargin{:});
