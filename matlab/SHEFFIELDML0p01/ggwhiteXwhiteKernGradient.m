function [g1, g2] = ggwhiteXwhiteKernGradient(ggwhiteKern, whiteKern, ...
    x, x2, covPar)

% GGWHITEXWHITEKERNGRADIENT Compute gradient between the GGWHITE and WHITE kernels.
%
%	Description:
%
%	[G1, G2] = GGWHITEXWHITEKERNGRADIENT(GGWHITEKERN, WHITEKERN, X)
%	computes the gradient of an objective function with respect to cross
%	kernel terms between GGWHITE and WHITE kernels for the multiple
%	output kernel.
%	 Returns:
%	  G1 - gradient of objective function with respect to kernel
%	   parameters of GGWHITE kernel.
%	  G2 - gradient of objective function with respect to kernel
%	   parameters of WHITE kernel.
%	 Arguments:
%	  GGWHITEKERN - the kernel structure associated with the GGWHITE
%	   kernel.
%	  WHITEKERN - the kernel structure associated with the WHITE kernel.
%	  X - inputs for which kernel is to be computed.
%
%	[G1, G2] = GGWHITEXWHITEKERNGRADIENT(GGWHITEKERN, WHITEKERN, X1, X2)
%	computes the gradient of an objective function with respect to cross
%	kernel terms between GGWHITE and WHITE kernels for the multiple
%	output kernel.
%	 Returns:
%	  G1 - gradient of objective function with respect to kernel
%	   parameters of GGWHITE kernel.
%	  G2 - gradient of objective function with respect to kernel
%	   parameters of WHITE kernel.
%	 Arguments:
%	  GGWHITEKERN - the kernel structure associated with the GGWHITE
%	   kernel.
%	  WHITEKERN - the kernel structure associated with the WHITE kernel.
%	  X1 - row inputs for which kernel is to be computed.
%	  X2 - column inputs for which kernel is to be computed.
%	whiteKernParamInit
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, GGWHITEKERNPARAMINIT, 


%	Copyright (c) 2008 Mauricio A. Alvarez and Neil D. Lawrence

   
if nargin < 5
    covPar = x2;
    x2 = x;
end

P = ggwhiteKern.precisionG;
Pinv = 1./P;
detPinv = prod(Pinv);
sqrtP = sqrt(P);
sqrtPx = x*sparseDiag(sqrtP);
sqrtPx2 = x2*sparseDiag(sqrtP);
n2 = dist2(sqrtPx, sqrtPx2);

factor = ggwhiteKern.variance*whiteKern.variance...
    /((pi*detPinv)^(ggwhiteKern.inputDimension/4));

factorSens = whiteKern.variance...
    /((pi*detPinv)^(ggwhiteKern.inputDimension/4));

factorVar = ggwhiteKern.variance...
    /((pi*detPinv)^(ggwhiteKern.inputDimension/4));

Kbase = exp(-0.5*n2);

k = factor*Kbase;


matGrad = zeros(ggwhiteKern.inputDimension,1);

for i = 1:ggwhiteKern.inputDimension,
    X = repmat(x(:,i),1, size(x2,1));
    X2 = repmat(x2(:,i)',size(x,1),1);
    if ggwhiteKern.isArd
        matGrad(i) = sum(sum(0.5*covPar.*k.*(Pinv(i) - (X - X2).*(X - X2))));
    else
        matGrad(i) = sum(sum(0.5*covPar.*k.*(Pinv - (X - X2).*(X - X2))));
    end
end

if ~ggwhiteKern.isArd
    matGrad = sum(matGrad);
end

g1 = [matGrad' 0 factorSens*sum(sum(covPar.*Kbase))];

g2 = factorVar*sum(sum(covPar.*Kbase));





