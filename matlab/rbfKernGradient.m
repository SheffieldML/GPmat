function g = rbfKernGradient(kern, x, varargin)

% RBFKERNGRADIENT Gradient of RBF kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% radial basis function
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
% SEEALSO rbfKernParamInit, kernGradient, rbfKernDiagGradient, kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006, 2009
%
% MODIFICATIONS : Mauricio Alvarez, 2009, David Luengo, 2009

% KERN


% The last argument is covGrad
if nargin < 4
  [k, sk, dist2xx] = rbfKernCompute(kern, x);
else
  [k, sk, dist2xx] = rbfKernCompute(kern, x, varargin{1});
end
% if gK is cell then return cell of gs
if isfield(kern, 'isNormalised') && (kern.isNormalised == true)
    g(1, 1) = (0.5*kern.variance/sqrt(2*pi)) * sum(sum(varargin{end} .* sk ...
        .* (1/sqrt(kern.inverseWidth)-sqrt(kern.inverseWidth)*dist2xx)));
    g(1, 2) = sqrt(kern.inverseWidth/(2*pi)) * sum(sum(varargin{end}.*sk));
else
    g(1, 1) = - 0.5 * sum(sum(varargin{end}.*k.*dist2xx));
    g(1, 2) = sum(sum(varargin{end}.*sk));
end
%/~
if any(isnan(g))
  warning('g is NaN')
end
%~/
