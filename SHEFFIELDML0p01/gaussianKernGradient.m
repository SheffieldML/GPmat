function g = gaussianKernGradient(kern, x, varargin)

% GAUSSIANKERNGRADIENT Gradient of gaussian kernel's parameters.
%
%	Description:
%
%	G = GAUSSIANKERNGRADIENT(KERN, X, PARTIAL) computes the gradient of
%	functions with respect to the gaussian kernel's parameters. As well
%	as the kernel structure and the input positions, the user provides a
%	matrix PARTIAL which gives the partial derivatives of the function
%	with respect to the relevant elements of the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters. The ordering of the vector should match that
%	   provided by the function kernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  X - the input locations for which the gradients are being
%	   computed.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The argument takes the
%	   form of a square matrix of dimension  numData, where numData is
%	   the number of rows in X.
%
%	G = GAUSSIANKERNGRADIENT(KERN, X1, X2, PARTIAL) computes the
%	derivatives as above, but input locations are now provided in two
%	matrices associated with rows and columns of the kernel matrix.
%	 Returns:
%	  G - gradients of the function of interest with respect to the
%	   kernel parameters.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are being
%	   computed.
%	  X1 - the input locations associated with the rows of the kernel
%	   matrix.
%	  X2 - the input locations associated with the columns of the kernel
%	   matrix.
%	  PARTIAL - matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The matrix should have
%	   the same number of rows as X1 and the same number of columns as X2
%	   has rows.
%	gaussianKernDiagGradient, kernGradX
%	
%	
%
%	See also
%	GAUSSIANKERNPARAMINIT, KERNGRADIENT, 


%	Copyright (c) 2008 Mauricio A. Alvarez and Neil D. Lawrence


%	With modifications by Mauricio A. Alvarez 2009.

  
if nargin < 4
    x2 = x;
    covPar = varargin{1};
else
    x2 = varargin{1};
    covPar = varargin{2};
end

if kern.isArd
    [K, Kbase] = gaussianKernCompute(kern, x, x2);
    matGrad = zeros(kern.inputDimension,1);
    temp = 0.5*covPar.*K;
    for i = 1:kern.inputDimension,
        pX = x(:,i);
        X = pX(:, ones(1, size(x2,1)));
        pX2 = x2(:,i)';
        X2 = pX2(ones(size(x,1),1), :);
%         X = repmat(x(:,i),1, size(x2,1));
%         X2 = repmat(x2(:,i)',size(x,1),1);
        X_X2 = X - X2;
        matGrad(i) = -sum(sum(temp.*X_X2.*X_X2));
    end
else
    [K, Kbase, n2] = gaussianKernCompute(kern, x, x2);   
    matGrad = - 0.5*sum(sum(covPar.*K.*n2,2));

end
g = [matGrad(:)' sum(sum(covPar.*Kbase))];





