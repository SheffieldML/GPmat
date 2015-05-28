function gX = ggwhiteXgaussianwhiteKernGradX(ggwhiteKern, gaussianwhiteKern, X, X2)

% GGWHITEXGAUSSIANWHITEKERNGRADX Compute gradient between the GG white and
% GAUSSIAN white kernels wrt the input locations
% FORMAT
% DESC computes the gradient between the 
%   GG white and GAUSSIAN white kernels with respect to the input positions where both
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
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009.

% KERN

if nargin < 3,    
    X2 = X;
else
    U = X;
    X = X2;
    X2 = U;
end

[K, P] = ggwhiteXgaussianwhiteKernCompute(ggwhiteKern, gaussianwhiteKern, X2, X);

if gaussianwhiteKern.nIndFunct == 1,    
    if ggwhiteKern.isArd
        PX = X*diag(P);
        PX2 = X2*diag(P);
    else
        PX = P*X;
        PX2 = P*X2;
    end
end

gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));

for i = 1:size(X, 1);
    if ggwhiteKern.isArd
        if gaussianwhiteKern.nIndFunct == 1,
            gX(:, :, i) = gaussianKernGradXpoint(PX(i, :), PX2, K(:,i));
        else
            for j = 1:size(X,2)
                partialDer = K(:,i).*P(:,i,j);
                gX(:, j, i) = gaussianKernGradXpoint( X(i, j), X2(:,j), partialDer);
            end
        end
    else
        if gaussianwhiteKern.nIndFunct == 1,
            gX(:, :, i) = gaussianKernGradXpoint(PX(i, :), PX2, K(:,i));
        else
            partialDer = K(:,i).*P(:,i);
            gX(:, :, i) = gaussianKernGradXpoint( X(i, :), X2, partialDer);
        end
    end
end

function gX = gaussianKernGradXpoint(x, x2, partialDer)

% GAUSSIANKERNGRADXPOINT Gradient with respect to one point of x.

gX = zeros(size(x2));
for i = 1:size(x, 2)
  gX(:, i) = (x2(:, i) - x(i)).*partialDer;
end


