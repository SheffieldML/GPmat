function gT = rbfinfwhiteKernGradX(kern, t1, t2)

% RBFINFWHITEKERNGRADX Gradient of RBF-WHITE kernel (with integration limits
% between minus infinity and infinity) with respect to a point t.
% FORMAT
% DESC computes the gradient of the RBF-WHITE kernel with respect to the
% input positions. 
% ARG kern : kernel structure for which gradients are being computed.
% ARG t1 : locations against which gradients are being computed.
% RETURN gT : the returned gradients. The gradients are returned in
% a matrix which is numData x numInputs x numData. Where numData is
% the number of data points and numInputs is the number of input
% dimensions in t1 (currently always one).
%
% FORMAT
% DESC computes the gradient of the RBF-WHITE kernel with respect to the
% input positions where both the row positions and column positions are
% provided separately.
% ARG kern : kernel structure for which gradients are being
% computed.
% ARG t1 : row locations against which gradients are being computed.
% ARG t2 : column locations against which gradients are being computed.
% RETURN gT : the returned gradients. The gradients are returned in
% a matrix which is numData2 x numInputs x numData1. Where numData1 is
% the number of data points in t1, numData2 is the number of data
% points in t2 and numInputs is the number of input dimensions in t1
% and t2 (currently always one).
%
% SEEALSO rbfinfwhiteKernParamInit, kernGradX, rbfinfwhiteKernDiagGradX
%
% COPYRIGHT : David Luengo, 2009

% KERN


if nargin < 3
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

gT = zeros(size(t1, 1), 1, size(t2, 1));

% Parameters of the kernel required in the computation
variance = kern.variance;
invWidth = kern.inverseWidth;

for i = size(t1, 1)
    gT(i, 1, :) = - ( variance * invWidth * (t1(i)-t2) * sqrt(invWidth/pi) / 4) ...
        .* exp(- invWidth * ((t1(i)-t2).^2) / 4);
end
