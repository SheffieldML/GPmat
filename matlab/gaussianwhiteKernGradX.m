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
% ARG x1 : row locations against which gradients are being computed.
% ARG x2 : column locations against which gradients are being computed.
%	
% SEEALSO : gaussianwhiteKernParamInit, kernGradX, gaussianwhiteKernDiagGradX
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008

% KERN

if nargin < 3,
    %covGrad = X2;
    X2 = X;
end

K  = gaussianwhiteKernCompute(kern, X, X2);
P = kern.precisionT/2;
PX = X*sparseDiag(P);
PX2 = X2*sparseDiag(P);

gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = gaussianwhiteKernGradXpoint(K(i,:)', PX(i, :), PX2);
end

%gXu = zeros(size(X));

if nargin <4,
    gX = gX*2;
    dgKX = gaussianwhiteKernDiagGradX(kern, X);
    for i = 1:size(X,1)
        gX(i, :, i) = dgKX(i, :);
    end
end

% for i = 1:size(X,1),
%     for j=1:size(X,2),
%       gXu(i,j) = covGrad(i,:)*gX(:,j,i);
%     end
% end

function gX = gaussianwhiteKernGradXpoint(gaussianPart, x, X2)

% GAUSSIANWHITEKERNGRADXPOINT Gradient with respect to one point of x.

gX = zeros(size(X2));
for i = 1:size(x, 2)
  gX(:, i) = (X2(:, i) - x(i)).*gaussianPart;
end
