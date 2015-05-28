function g = dexpKernGradient(kern, x1, varargin)

% DEXPKERNGRADIENT Gradient of the double exponential kernel's parameters.
%
% FORMAT
% DESC computes the gradient of functions with respect to the double
% exponential kernel's parameters. As well as the kernel structure
% and the input positions, the user provides a matrix PARTIAL which gives
% the partial derivatives of the function with respect to the relevant
% elements of the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x1 : the input locations for which the gradients are being
% computed in the form of a design matrix.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in x1.
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
% kernel matrix in the form of a design matrix.
% ARG x2 : the input locations associated with the columns of the
% kernel matrix in the form of a design matrix.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as x1 and the same number of columns
% as x2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
%
% SEEALSO dexpKernParamInit, kernGradient, dexpKernDiagGradient, kernGradX
%
% COPYRIGHT : David Luengo, 2009

% KERN


if nargin < 4
    x2 = x1;
else
    x2 = varargin{1};
end

[K, sK, n1] = dexpKernCompute(kern, x1, x2);

g(1) = kern.variance * sum(sum(varargin{end} .* sK ...
    .* (1 - n1*kern.decay)));
g(2) = kern.decay * sum(sum(varargin{end} .* sK));

%/~
if any(isnan(g))
  warning('g is NaN')
end
%~/
