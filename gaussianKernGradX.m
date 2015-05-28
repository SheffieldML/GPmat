function gX = gaussianKernGradX(kern,X, X2, covGrad)

% GAUSSIANKERNGRADX Gradient of gaussian kernel with respect to input locations.
% FORMAT
% DESC computes the gradient of the
%	gaussian kernel with respect to the input positions where both
%	the row positions and column positions are provided separately.
% RETURN g : the returned gradients. 
% ARG kern : kernel structure for which gradients are being computed.
% ARG x1 : row locations against which gradients are being computed.
% ARG x2 : column locations against which gradients are being computed.
% ARG covGrad : partial derivatives
%	
% SEEALSO : gaussianKernParamInit, kernGradX, gaussianKernDiagGradX
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009, 2010

% KERN


if nargin < 4,
    covGrad = X2;
    X2 = X;
end
K = gaussianKernCompute(kern, X, X2);
P = kern.precisionU;
if kern.isArd
    PX = X*diag(P);
    PX2 = X2*diag(P);
else
    PX = P*X;
    PX2 = P*X2;
end
gX = zeros(size(X));

temp = (covGrad.*K)';

for i=1:size(X,2),
    mPX2 = PX2(:,i);
    MPX2 = mPX2(:,ones(1,size(X,1)));
    mPX = PX(:,i)';
    MPX = mPX(ones(size(X2,1),1),:);
    gX(:,i) = sum(temp.*(MPX2-MPX),1)';
end

if nargin <4,
    gX = gX*2;
end

%/~ Old version of the code
% if nargin < 3,
%     X2 = X;
% end
% 
% K = gaussianKernCompute(kern, X, X2);
% 
% P = kern.precisionU;
% 
% if kern.isArd
%     PX = X*diag(P);
%     PX2 = X2*diag(P);
% else
%     PX = P*X;
%     PX2 = P*X2;
% end
% 
% gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
% for i = 1:size(X, 1);
%   gX(:, :, i) = gaussianKernGradXpoint(K(i,:)', PX(i, :), PX2);
% end
% 
% if nargin <3,
%     gX = gX*2;
%     dgKX = gaussianKernDiagGradX(kern, X);
%     for i = 1:size(X,1)
%         gX(i, :, i) = dgKX(i, :);
%     end
% end
% 
% function gX = gaussianKernGradXpoint(gaussianPart, x, X2)
% 
% % RBFKERNGRADXPOINT Gradient with respect to one point of x.
% 
% gX = zeros(size(X2));
% for i = 1:size(x, 2)
%   gX(:, i) = (X2(:, i) - x(i)).*gaussianPart;
% end
%~/
