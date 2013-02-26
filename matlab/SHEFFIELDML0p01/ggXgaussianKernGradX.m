function gX = ggXgaussianKernGradX(ggKern, gaussianKern, X, X2, covGrad)

% GGXGAUSSIANKERNGRADX Compute gradient between the GG and GAUSSIAN
%
%	Description:
%	kernels wrt the input locations
%
%	G = GGXGAUSSIANKERNGRADX(KERN, X1, X2) computes the gradient between
%	the GG and GAUSSIAN kernels with respect to the input positions
%	where both the row positions and column positions are provided
%	separately.
%	 Returns:
%	  G - the returned gradients.
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  X1 - row locations against which gradients are being computed.
%	  X2 - column locations against which gradients are being computed.
%	
%	
%
%	See also
%	GAUSSIANKERNPARAMINIT, KERNGRADX, GAUSSIANKERNDIAGGRADX


%	Copyright (c) 2008 Mauricio A. Alvarez and Neil D. Lawrence


%	With modifications by Mauricio A. Alvarez 2009, 2010


if nargin < 5,
    covGrad = X2;
    X2 = X;
else
    U = X;
    X = X2;
    X2 = U;
end


[K, Kbase, Pqrinv, Prinv, P] = ggXgaussianKernCompute(ggKern, ....
    gaussianKern, X2, X);

if ggKern.isArd
    PX = X*diag(P);
    PX2 = X2*diag(P);
else
    PX = P*X;
    PX2 = P*X2;
end

temp = covGrad'.*K;

gX2 = zeros(size(X));
for i=1:size(X,2),
    mPX2 = PX2(:,i);
    MPX2 = mPX2(:,ones(1,size(X,1)));
    mPX = PX(:,i)';
    MPX = mPX(ones(size(X2,1),1),:);
    gX2(:,i) = sum(temp.*(MPX2-MPX),1)';
end

gX = gX2;


