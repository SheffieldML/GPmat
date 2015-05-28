function gX = dexpKernGradX(kern, x1, x2)

% DEXPKERNGRADX Gradient of the double exponential kernel with respect to a
% point x.
%
% FORMAT
% DESC computes the gradient of the double exponential kernel kernel with
% respect to the input positions.
% ARG kern : kernel structure for which gradients are being computed.
% ARG x1 : locations against which gradients are being computed in the form
% of a design matrix.
% RETURN gX : the returned gradients. The gradients are returned in
% a matrix which is numData x numInputs x numData. Where numData is
% the number of data points and numInputs is the number of input
% dimensions in x1.
%
% FORMAT
% DESC computes the gradident of the double exponential kernel with respect
% to the input positions where both the row positions and column positions
% are provided separately.
% ARG kern : kernel structure for which gradients are being computed.
% ARG x1 : row locations against which gradients are being computed in the
% form of a design matrix.
% ARG x2 : column locations against which gradients are being computed in
% the form of a design matrix.
% RETURN gX : the returned gradients. The gradients are returned in
% a matrix which is numData2 x numInputs x numData1. Where numData1 is
% the number of data points in x1, numData2 is the number of data
% points in x2 and numInputs is the number of input dimensions in x1
% (currently always one).
%
% SEEALSO dexpKernParamInit, kernGradX, dexpKernDiagGradX
%
% COPYRIGHT : David Luengo, 2009

% KERN


if nargin < 3
  x2 = x1;
end
[Nx1, dx1] = size(x1);
[Nx2, dx2] = size(x2);
if (dx1 ~= dx2)
    error('The dimension of x1 and x2 must be the same');
end

% Initialisation of the gradient matrix
gX = zeros(Nx1, dx1, Nx2);

% Computation
K = dexpKernCompute(kern, x1, x2);
for i = size(x1, 1)
    gX(i, :, :) = - kern.decay * transpose(sign(repmat(x1(i, :), Nx2, 1)-x2)) ...
        .* repmat(K(i, :), dx1, 1);
end
