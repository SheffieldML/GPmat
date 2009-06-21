function gX = gaussianwhiteKernGradX(kern,X, X2)

% GAUSSIANWHITEKERNGRADX Gradient of gaussian white kernel with respect 
%                        to input locations.
% FORMAT
% DESC computes the gradient of the
%	gaussian white kernel with respect to the input positions where both
%	the row positions and column positions are provided separately.
% RETURN g : the returned gradients. The gradients are returned in a matrix
%	   which is numData2 x numInputs x numData1. Where numData1 is the
%	   number of data points in X1, numData2 is the number of data points
%	   in X2 and numInputs is the number of input dimensions in X.
% ARG kern : kernel structure for which gradients are being computed.
% ARG X  : row locations against which gradients are being computed.
% ARG X2 : column locations against which gradients are being computed.
%	
% SEEALSO : gaussianwhiteKernParamInit, kernGradX, gaussianwhiteKernDiagGradX
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009.

% KERN

if nargin < 3,
    X2 = X;
end

if kern.isArd
    K  = gaussianwhiteKernCompute(kern, X, X2);
    P = kern.precisionT/2;
    PX = X*sparseDiag(P);
    PX2 = X2*sparseDiag(P);
    gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
    for i = 1:size(X, 1);
        gX(:, :, i) = gaussianwhiteKernGradXpoint(PX(i, :), PX2, K(i,:)');
    end
else
    precCols = repmat(kern.precisionT', 1, size(X,1));
    precRows = repmat(kern.precisionT , size(X2,1), 1);
    precColsInv = 1./precCols;
    precRowsInv = 1./precRows;
    Pinv = precColsInv + precRowsInv;
    P = 1./Pinv;
    detPinv = Pinv.^kern.inputDimension;
    factor = 1./((2*pi)^(kern.inputDimension/2)*sqrt(detPinv));
    K = kern.sigma2Noise*factor.*exp(-0.5*P.*dist2(X, X2));
    gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
    for i = 1:size(X, 1);
        partialDer = (K(i,:).*P(i,:))';
        gX(:, :, i) = gaussianwhiteKernGradXpoint( X(i, :), X2, partialDer);
    end
end
gX = gX*2;
dgKX = gaussianwhiteKernDiagGradX(kern, X);
for i = 1:size(X,1)
    gX(i, :, i) = dgKX(i, :);
end

function gX = gaussianwhiteKernGradXpoint(x, x2, partialDer)

% GAUSSIANWHITEKERNGRADXPOINT Gradient with respect to one point of x.

gX = zeros(size(x2));
for i = 1:size(x, 2)
  gX(:, i) = (x2(:, i) - x(i)).*partialDer;
end
