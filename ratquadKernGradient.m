function g = ratquadKernGradient(kern, x, varargin)

% RATQUADKERNGRADIENT Gradient of RATQUAD kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% rational quadratic
% kernel's parameters. As well as the kernel structure and the
% input positions, the user provides a matrix PARTIAL which gives
% the partial derivatives of the function with respect to the
% relevant elements of the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x : the input locations for which the gradients are being
% computed. 
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in X.
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
% ARG x1 : the input locations associated with the rows of the
% kernel matrix.
% ARG x2 : the input locations associated with the columns of the
% kernel matrix.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as X1 and the same number of columns
% as X2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
%
% SEEALSO ratquadKernParamInit, kernGradient, ratquadKernDiagGradient, kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN


% The last argument is covGrad
wi2 = .5/(kern.lengthScale*kern.lengthScale*kern.alpha);
if nargin < 4
  n2 = dist2(x, x);
else
  n2 = dist2(x, varargin{1});
end
baseVal = (1+n2*wi2);
kbase = baseVal.^-kern.alpha;
kbase2 = -kern.alpha*kbase./baseVal;
g(1) = -kern.variance*sum(sum(varargin{end}.*(n2*wi2/kern.alpha.*kbase2 ...
                                             + log(baseVal).*kbase)));
g(2) = -kern.variance*sum(sum(varargin{end}.*n2.*kbase2))*2*wi2/kern.lengthScale;
g(3) = sum(sum(varargin{end}.*(1+n2*wi2).^-kern.alpha));
