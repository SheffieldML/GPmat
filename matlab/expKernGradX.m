function gX = expKernGradX(kern, x, x2)

% EXPKERNGRADX Gradient of EXP kernel with respect to a point x.
% FORMAT
% DESC computes the gradient of the exponentiated
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
% DESC computes the gradident of the exponentiated
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
% SEEALSO expKernParamInit, kernGradX, expKernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

[K, argK, kDiagdMat] = expKernCompute(kern, x, x2);
gX = kernGradX(kern.argument, x, x2);
if kern.isStationary
  for i = 1:size(gX, 3)
    gX(:, :, i) = gX(:, :, i)...
        .*repmat(exp(argK(i, :)'), 1, size(x, 2))...
        .*exp(kDiagMat)*kern.variance;
  end
else

%   for i = 1:size(gX, 3)
%     gX(:, :, i) = gX(:, :, i) ...
%         .*(repmat(exp(argK(i, :)').*exp(kDiagMat(i, :)'), 1, size(x, 2))...
%         + repmat(((exp(argK(i, :)')-1).*exp(kDiagMat(i, :))*kern.variance;
%   end
    
end
