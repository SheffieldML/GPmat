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
    sqrtP = sqrt(kern.precisionT/2);
    Px = x*sparseDiag(sqrtP);   
    if nargin < 4
        n2 = dist2(Px, Px);
    else
        Px2 = x2*sparseDiag(sqrtP);
        n2 = dist2(Px, Px2);
    end
    kBase = exp(-0.5*n2);
    k = kern.sigma2Noise*kBase;
    matGrad = zeros(kern.inputDimension,1);
    for i = 1:kern.inputDimension,
        X = repmat(x(:,i),1, size(x2,1));
        X2 = repmat(x2(:,i)',size(x,1),1);
        matGrad(i) = sum(sum(0.25*covPar.*k.*(- (X - X2).*(X - X2))));
    end
    g = [matGrad' sum(sum(covPar.*kBase))];
else
    dist = dist2(x, x2);
    if nargin < 4
        if kern.nIndFunct~=size(x,1)
            error(['The number of inducing functions must be equal the' ...
                'number of inducing points']);
        end
        precCols = repmat(kern.precisionT', 1, size(x,1));
        precRows = repmat(kern.precisionT , size(x,1), 1);
        precColsInv = 1./precCols;
        precRowsInv = 1./precRows;
        Pinv = precColsInv + precRowsInv;
        P = 1./Pinv;
        precColInvNum =  repmat((kern.precisionT.^(-kern.inputDimension/4))', 1, size(x,1));
        precRowsInvNum=  repmat(kern.precisionT.^(-kern.inputDimension/4) , size(x,1), 1);
        factorDen = Pinv.^(kern.inputDimension/2); 
        factor = 2^(kern.inputDimension/2)*(precColInvNum.*precRowsInvNum)./factorDen;
        n2  = P.*dist;
        kBase = exp(-0.5*n2);
        k = kern.sigma2Noise*factor.*kBase;
        gradSigma2Noise = sum(sum(factor.*covPar.*kBase));
    else
        precCols = kern.precisionT'/2;
        precColsInv = 1./precCols;
        Pinv =  precColsInv;
     %   P = repmat(1./Pinv', size(x,1), 1);    
        P = repmat(1./Pinv, 1, size(x2,1));
        n2  = P.*dist;
        kBase = exp(-0.5*n2);
        k = kern.sigma2Noise*kBase;
        gradSigma2Noise = sum(sum(covPar.*kBase));
    end   
    if nargin<4
        temp = ((P.*precColsInv).^2).*(kern.inputDimension*Pinv - dist);
        matGrad = 2*sum(0.5*covPar.*k.*(-0.5*kern.inputDimension*precColsInv+temp), 2)';
    else
        matGrad = 0.25*sum(covPar.*k.*( - dist),2)';
    end
    g = [matGrad gradSigma2Noise];
end



