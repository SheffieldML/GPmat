function g = gaussianwhiteKernGradient(kern, x, varargin)

% GAUSSIANWHITEKERNGRADIENT Gradient of gaussian white kernel's parameters.
% FORMAT
% DESC computes the gradient of
%	functions with respect to the gaussian white kernel's
%	parameters. As well as the kernel structure and the input positions,
%	the user provides a matrix PARTIAL which gives the partial
%	derivatives of the function with respect to the relevant elements of
%	the kernel matrix.
% RETURN g:  gradients of the function of interest with respect to the
%	   kernel parameters. The ordering of the vector should match that
%	   provided by the function kernExtractParam.
% ARG kern : the kernel structure for which the gradients are being
%	   computed.
% ARG x : the input locations for which the gradients are being
%	   computed.
% ARG partial : matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The argument takes the
%	   form of a square matrix of dimension  numData, where numData is
%	   the number of rows in X.
%
% FORMAT
% DESC  computes the derivatives
%	as above, but input locations are now provided in two matrices
%	associated with rows and columns of the kernel matrix.
% RETURN g : gradients of the function of interest with respect to the
%	   kernel parameters.
% ARG kern : the kernel structure for which the gradients are being
%	   computed.
% ARG x1 : the input locations associated with the rows of the kernel
%	   matrix.
% ARG x2 : the input locations associated with the columns of the kernel
%	   matrix.
% ARG partial : matrix of partial derivatives of the function of
%	   interest with respect to the kernel matrix. The matrix should have
%	   the same number of rows as X1 and the same number of columns as X2
%	   has rows.
% 
% SEEALSO : gaussianwhiteKernParamInit, kernGradient,
%           gaussianwhiteKernDiagGradient, kernGradX
%  
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009

% KERN
  
if nargin < 4
    x2 = x;
    covPar = varargin{1};
else
    x2 = varargin{1};
    covPar = varargin{2};
end

if kern.isArd
    if kern.nIndFunct == 1
        [k, P, kBase] = gaussianwhiteKernCompute(kern, x, x2);
        matGrad = zeros(kern.inputDimension,1);
        for i = 1:kern.inputDimension,
            X = repmat(x(:,i),1, size(x2,1));
            X2 = repmat(x2(:,i)',size(x,1),1);
            matGrad(i) = sum(sum(0.25*covPar.*k.*(- (X - X2).*(X - X2))));
        end
        gradSigma2Noise = sum(sum(covPar.*kBase));       
    else
        matGrad = zeros(kern.inputDimension, kern.nIndFunct);
        if nargin < 4,
            [K, P, Kbase, dist, Factor, precColsInv, Pinv] = ...
                gaussianwhiteKernCompute(kern, x);
             gradSigma2Noise = sum(sum(Factor.*covPar.*Kbase)); 
        else
            [K, P, Kbase, dist]  = gaussianwhiteKernCompute(kern, x, x2);
            gradSigma2Noise = sum(sum(covPar.*Kbase));
        end        
        for i=1:kern.inputDimension,
            if nargin<4
                temp = ((P(:,:,i).*precColsInv(:,:,i)).^2).*(Pinv(:,:,i) - dist(:,:,i));
                matGrad(i,:) = 2*sum(0.5*covPar.*K.*(-0.5*precColsInv(:,:,i)+temp), 2)';
            else
                matGrad(i,:) = 0.25*sum(covPar.*K.*( - dist(:,:,i)),2)';
            end            
        end
    end
else
    if kern.nIndFunct == 1
        if nargin<4
            x2 = x;
        end
        [K, P, Kbase, n2] = gaussianwhiteKernCompute(kern, x, x2);
        matGrad = 0.25*sum(sum(covPar.*K.*( - n2),2));
        gradSigma2Noise = sum(sum(covPar.*Kbase));
    else        
        if nargin<4
            [K, P, Kbase, dist, Factor, precColsInv, Pinv] = ...
                gaussianwhiteKernCompute(kern, x);
            temp = ((P.*precColsInv).^2).*(kern.inputDimension*Pinv - dist);
            matGrad = 2*sum(0.5*covPar.*K.*(-0.5*kern.inputDimension*precColsInv+temp), 2)';
            gradSigma2Noise = sum(sum(Factor.*covPar.*Kbase));
        else
            [K, P, Kbase, dist]  = gaussianwhiteKernCompute(kern, x, x2);
            matGrad = 0.25*sum(covPar.*K.*( - dist),2)';
            gradSigma2Noise = sum(sum(covPar.*Kbase));
        end
    end
end
g = [matGrad(:)' gradSigma2Noise];
