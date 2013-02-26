function gX = gaussianKernGradX(kern,X, X2, covGrad)

% GAUSSIANKERNGRADX Gradient of gaussian kernel with respect to input locations.
%
%	Description:
%
%	G = GAUSSIANKERNGRADX(KERN, X1, X2, COVGRAD) computes the gradient
%	of the gaussian kernel with respect to the input positions where
%	both the row positions and column positions are provided separately.
%	 Returns:
%	  G - the returned gradients.
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  X1 - row locations against which gradients are being computed.
%	  X2 - column locations against which gradients are being computed.
%	  COVGRAD - partial derivatives
%	
%	
%
%	See also
%	GAUSSIANKERNPARAMINIT, KERNGRADX, GAUSSIANKERNDIAGGRADX


%	Copyright (c) 2008 Mauricio A. Alvarez and Neil D. Lawrence


%	With modifications by Mauricio A. Alvarez 2009, 2010



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

