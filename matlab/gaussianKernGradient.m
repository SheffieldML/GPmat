function g = gaussianKernGradient(kern, x, varargin)

% GAUSSIANKERNGRADIENT Gradient of gaussian kernel's parameters.
% FORMAT
% DESC computes the gradient of
%	functions with respect to the gaussian kernel's
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
% SEEALSO : gaussianKernParamInit, kernGradient,
% gaussianKernDiagGradient, kernGradX
%  
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008

% KERN
  
if nargin < 4
    x2 = x;
    covPar = varargin{1};
else
    x2 = varargin{1};
    covPar = varargin{2};
end

L = sqrt(kern.precision_u);
Lx = x*diag(L);
Lx2 = x2*diag(L);
n2 = dist2(Lx, Lx2);


if isfield(kern, 'isNormalised') && ~isempty(kern.isNormalised)
    if kern.isNormalised
        option = 1;
    else
        option = 0;
    end
else
   option = 0; 
end

if option
    kBase = sqrt(prod(kern.precision_u))*exp(-0.5*n2);
    k = kern.sigma2_u*kBase;
    matGrad = zeros(kern.inputDimension,1);

    for i = 1:kern.inputDimension,
        X = repmat(x(:,i),1, size(x2,1));
        X2 = repmat(x2(:,i)',size(x,1),1);
        matGrad(i) = sum(sum(0.5*covPar.*k.*(1/(kern.precision_u(i)) - (X - X2).*(X - X2))));
    end
else
    kBase = exp(-0.5*n2);
    k = kern.sigma2_u*kBase;
    matGrad = zeros(kern.inputDimension,1);
    for i = 1:kern.inputDimension,
        X = repmat(x(:,i),1, size(x2,1));
        X2 = repmat(x2(:,i)',size(x,1),1);
        matGrad(i) = -sum(sum(0.5*covPar.*k.*(X - X2).*(X - X2)));
    end
end

g = [matGrad' sum(sum(covPar.*kBase))];





