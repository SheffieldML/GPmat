function g = matern52KernGradient(kern, x, varargin)

% MATERN52KERNGRADIENT Gradient of MATERN52 kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% matern kernel with nu=5/2
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
% SEEALSO matern52KernParamInit, kernGradient, matern52KernDiagGradient, kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2009

% KERN

% The last argument is covGrad
if nargin < 4
  n2 = dist2(x, x);
else
  n2 = dist2(x, varargin{1});
end
warnState = warning('query', 'MATLAB:divideByZero');
warning('off', 'MATLAB:divideByZero');
wi2 = (5/(kern.lengthScale*kern.lengthScale));
n2wi2 = n2*wi2;
sqrtn2wi2 = sqrt(n2wi2);
sk = (1+sqrtn2wi2+n2wi2/3).*exp(-sqrtn2wi2);
K = kern.variance*sk;
ratio = n2./sqrtn2wi2;
ratio(find(isnan(ratio)))=1;
g(1) = sum(sum(varargin{end}.*(wi2/kern.lengthScale*ratio.*(K-kern.variance.*exp(-sqrtn2wi2))-2*wi2/kern.lengthScale*kern.variance/3*n2.*exp(-sqrtn2wi2))));
warning(warnState.state, 'MATLAB:divideByZero');
g(2) =  sum(sum(varargin{end}.*sk));
