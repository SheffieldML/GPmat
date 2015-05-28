function g = ndsimKernGradient(kern, t1, varargin)

% NDSIMKERNGRADIENT Gradient of SIM kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% single input motif
% kernel's parameters. As well as the kernel structure and the
% input positions, the user provides a matrix PARTIAL which gives
% the partial derivatives of the function with respect to the
% relevant elements of the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG t : the input locations for which the gradients are being
% computed. 
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in X.
%
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters. The ordering of the vector should match
% that provided by the function kernExtractParam.
%
% FORMAT
% DESC computes the derivatives as above, but input locations are
% now provided in two matrices associated with rows and columns of
% the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG t1 : the input locations associated with the rows of the
% kernel matrix.
% ARG t2 : the input locations associated with the columns of the
% kernel matrix.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as X1 and the same number of columns
% as X2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
%
% That is, this function computes sum(sum(deriv(K,a).*partial)) for
% each scalar parameter a, where deriv(K,a) is the matrix-valued
% derivative of the kernel K (computed between times t and t2) with
% respect to the parameter a.
%
% SEEALSO simKernParamInit, kernGradient, simKernDiagGradient, kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Jaakko Peltonen, 2011

% KERN

if length(varargin)<2
  t2 = t1;
  covGrad = varargin{1};
else
  t2 = varargin{1};
  covGrad = varargin{2};
end



sigma = sqrt(2/kern.inverseWidth);
if isfield(kern, 'isNegativeS') && (kern.isNegativeS == true)
    variancemultiplier = (kern.sensitivity*kern.sensitivity);
else
    variancemultiplier = kern.variance;
end

dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);

k1a = sqrt(pi)/2*( t1Mat.*erf(t1Mat/sigma) + t2Mat.*erf(t2Mat/sigma) - diffT.*erf(diffT/sigma) );
k1b = sigma/2*( exp(-(t1Mat/sigma).^2) + exp(-(t2Mat/sigma).^2) - exp(-(diffT/sigma).^2) - 1 );
K = (k1a+k1b)*sigma*variancemultiplier;


% k = ndsimKernCompute(kern, t1, t2);
% deriv_variancemultiplier = sum(sum((k/variancemultiplier).*covGrad));
% deriv_sigma = k/sigma + variancemultiplier*sigma/2*(exp(-(t1Mat/sigma).^2)+exp(-(t2Mat/sigma).^2)-exp(-(diffT/sigma).^2)-1);
% deriv_sigma = sum(sum(deriv_sigma.*covGrad));
% % sigma=sqrt(2/inversewidth) 
% % --> df/dinversewidth = (df/dsigma)*(dsigma/dinversewidth)
% % = (df/dsigma)*((-1/2)*sqrt(2)*(inversewidth^(-3/2)))
% % = (df/dsigma)*(-1/2*sigma/inversewidth)
% deriv_inversewidth=deriv_sigma*(-0.5*sigma/kern.inverseWidth);


deriv_variancemultiplier = sum(sum(sigma*(k1a+k1b).*covGrad));
deriv_sigma = (k1a+2*k1b)*variancemultiplier;
% deriv_sigma = (k1a+k1b)*variancemultiplier ...
%     + variancemultiplier*sigma/2*(exp(-(t1Mat/sigma).^2)+exp(-(t2Mat/sigma).^2)-exp(-(diffT/sigma).^2)-1) ...
%     + variancemultiplier/sigma*(-(t1Mat.^2).*exp(-(t1Mat/sigma).^2)-(t2Mat.^2).*exp(-(t2Mat/sigma).^2)+(diffT.^2).*exp(-(diffT/sigma).^2)) ...
%    + variancemultiplier/sigma*(-(t1Mat.^2).*exp(-(t1Mat/sigma).^2)-(t2Mat.^2).*exp(-(t2Mat/sigma).^2)+(diffT.^2).*exp(-(diffT/sigma).^2)) ;
    
deriv_sigma = sum(sum(deriv_sigma.*covGrad));
deriv_inversewidth=deriv_sigma*(-0.5*sigma/kern.inverseWidth);


if isfield(kern, 'isNegativeS') && (kern.isNegativeS == true)
  % variancemultiplier=sensitivity^2
  % --> df/dsensitivity = (df/dvarmult)*(dvarmult/dsensitivity)
  % = (df/dvarmult)*(2*sensitivity)  
  deriv_sensitivity=deriv_variancemultiplier*2*kern.sensitivity;  
  g = [deriv_inversewidth deriv_sensitivity];
else
  g = [deriv_inversewidth deriv_variancemultiplier];
end


% gaussianinitial currently unsupported

