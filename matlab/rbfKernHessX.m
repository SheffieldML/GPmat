function hessX = rbfKernHessX(kern, X, X2)

% RBFKERNHESSX : Hessian Matrix of RBF kernel with respect to input 
% locations.
%
% FORMAT
% DESC computes the second derivatives of the radial basis function
% kernel with respect to the input positions where both the row
% positions and column positions are provided separately.
% ARG kern : kernel structure for which second derivatives are being
% computed.
% ARG x1 : row locations against which second derivatives are being 
% computed.
% ARG x2 : column locations against which second derivatives are being 
% computed.
% RETURN hessX : the returned Hessians. The Hessians are returned in
% a matrix which is numData2 x numInputs x numImputs x numData1, 
% where numData1 is the number of data points in X, numData2 is the number
% of data points in X2 and numInputs is the number of input dimensions in X.
%
% SEEALSO : rbfKernParamInit, kernGradX, rbfKernGradX, rbfKernDiagGradX
%
% COPYRIGHT : Neil Lawrence, 2004, 2005, 2006
%             Alessandra Tosi, 2013
%
% SHEFFIELDML

hessX = zeros(size(X2, 1), size(X2, 2), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  hessX(:, :, :,i) = rbfKernHessXpoint(kern, X(i, :), X2);
end
  

function hX = rbfKernHessXpoint(kern, x, X2)

% RBFKERNHESSXPOINT : Hessian Matrix with respect to one point x.

hX = zeros([size(X2, 1), size(X2, 2), size(X2, 2)]);
n2 = dist2(X2, x);
wi2 = (0.5 * kern.inverseWidth);
rbfPart = kern.variance*exp(-n2*wi2);
for i = 1:size(x, 2)
    for j = 1:size(x,2)
        if i==j
            hX(:,i,j) = - kern.inverseWidth*(kern.inverseWidth*(X2(:, i) - x(j)).^2 - 1).*rbfPart;
        else
            hX(:,i,j) = - kern.inverseWidth^2*(X2(:, i) - x(i)).*(X2(:, j) - x(j)).*rbfPart;
        end
    end
end
if isfield(kern, 'isNormalised') && (kern.isNormalised == true)
    hX = hX * sqrt(kern.inverseWidth/(2*pi));
end
