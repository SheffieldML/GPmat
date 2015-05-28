function g = ggKernGradient(kern, x, varargin)

% GGKERNGRADIENT Gradient of GG kernel's parameters.
% FORMAT
% DESC computes the gradient of
%	functions with respect to the gaussian gaussian kernel's
%	parameters. As well as the kernel structure and the input positions,
%	the user provides a matrix PARTIAL which gives the partial
%	derivatives of the function with respect to the relevant elements of
%	the kernel matrix.
% RETURN g : gradients of the function of interest with respect to the
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
% DESC computes the derivatives
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
% SEEALSO : ggKernParamInit, kernGradient, ggKernDiagGradient, kernGradX
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008

% KERN

if length(varargin)<2
  x2 = x;
else
  x2 = varargin{1};
end

covGrad = varargin{end};

if kern.isArd
    [K, Kbase, Prinv, Pqrinv, P] = ggKernCompute(kern, x, x2);
    matGradPqr = zeros(kern.inputDimension,1);
    matGradPr  = zeros(kern.inputDimension,1);
    temp = covGrad.*K;
    for i=1:kern.inputDimension
        pX = x(:,i);
        X = pX(:, ones(1, size(x2,1)));
        pX2 = x2(:,i)';
        X2 = pX2(ones(size(x,1),1), :);
        %X = repmat(x(:,i),1, size(x2,1));
        %X2 = repmat(x2(:,i)',size(x,1),1);
        X_X2 = (X - X2).*(X - X2);
        matGradPqr(i) = -sum(sum(temp.*(Pqrinv(i)*P(i)*X_X2*P(i)*Pqrinv(i))));
        matGradPr(i) = -0.5*sum(sum(temp.*(Prinv(i)*P(i)*X_X2*P(i)*Prinv(i))));
    end
else
    [K, Kbase, Prinv, Pqrinv, P, dist] = ggKernCompute(kern, x, x2);
    matGradPqr = -sum(sum(covGrad.*K.*((Pqrinv*P)^2*dist)));
    matGradPr  = -0.5*sum(sum(covGrad.*K.*((Prinv*P)^2*dist)));
end

gradSigma2Latent =  kern.sensitivity^2*sum(sum(covGrad.*Kbase));
gradSensitivity  = 2*kern.sensitivity*kern.sigma2Latent*sum(sum(covGrad.*Kbase)); 

g = [matGradPr(:)' matGradPqr(:)' gradSigma2Latent gradSensitivity];
