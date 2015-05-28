function g = kernGradX(kern, x, x2)

% KERNGRADX Compute the gradient of the kernel wrt X.
% FORMAT
% DESC computes the gradient of the kernel with respect to the
% input positions. 
% ARG kern : kernel structure for which gradients are being
% computed.
% ARG x : locations against which gradients are being computed.
% RETURN g : the returned gradients. The gradients are returned in
% a matrix which is numData x numInputs x numData. Where numData is
% the number of data points and numInputs is the number of input
% dimensions in X.
%
% FORMAT
% DESC computes the gradient of the kernel with respect to the
% input positions where both the row positions and column positions
% are provided separately.
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
% SEEALSO : kernDiagGradX, kernGradient
%
% COPYRIGHT : Neil D.Lawrence, 2004, 2005, 2006

% KERN

fhandle = str2func([kern.type 'KernGradX']);
if nargin < 3
  g = kernGradX(kern, x, x);
  dg = kernDiagGradX(kern, x);
  for i = 1:size(x, 1)
    g(i, :, i) = dg(i, :);
  end
else
  g = fhandle(kern, x, x2);
end
