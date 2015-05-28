function g = gibbsKernGradient(kern, x, varargin)

% GIBBSKERNGRADIENT Gradient of GIBBS kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% Mark Gibbs's non-stationary
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
% SEEALSO gibbsKernParamInit, kernGradient, gibbsKernDiagGradient, kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2009

% KERN

%  pretty(simple(diff((2*(l_i*l_j)/(l_i*l_i + l_j*l_j))^(d/2)*exp(-r*r/(l_i*l_i + l_j*l_j)),l_i)))
 
%                                         2
%       /    l_i l_j  \(1/2 d)           r               4        4        2  2
%   1/2 |2 -----------|        exp(- -----------) (-d l_i  + d l_j  + 4 l_i  r )
%       |     2      2|                 2      2
%       \  l_i  + l_j /              l_i  + l_j

%            /      2      2 2
%           /  ((l_i  + l_j )  l_i)
%          /
% >> pretty(simple(diff((2*(l_i*l_i)/(l_i*l_i + l_i*l_i))^(d/2)*exp(-r*r/(l_i*l_i + l_i*l_i)),l_i)))
 
%                                              2
%                                2            r
%                               r  exp(- 1/2 ----)
%                                               2
%                                            l_i
%                               ------------------
%                                         3
%                                      l_i
% >> pretty(simple(diff((2*(l_i*l_i)/(l_i*l_i + l_i*l_i))^(d/2)*exp(-r*r/(l_i*l_i + l_i*l_i)),l_i)))

fhandle = str2func([kern.lengthScaleTransform, 'Transform']);
g = zeros(1, kern.nParams);
% The last argument is covGrad
if length(varargin)<2
  [k, sk, n2, w2, l] = gibbsKernCompute(kern, x);
  gOut = modelOutputGrad(kern.lengthScaleFunc, x);
  gradFact = fhandle(l, 'gradfact');
  L1 = repmat(l, 1, size(l, 1));
  L2 = L1';
  covGrad = varargin{end};
  covGrad(1:size(covGrad, 1)+1:end) = 0;
  base = covGrad.*k;
  base2 = base.*(kern.inputDimension/2*(L2.^4 - L1.^4) + 2*L1.*L1.*n2)./(w2.*w2.*L1);
  for i = 1:size(g, 2)-1
    g(i) = g(i) + 2*sum(sum(base2.*repmat(gOut(:, i).*gradFact, 1, size(x, 1))));
  end
%/~
%  covGrad = varargin{end};
%  covGradDiag = diag(covGrad);
%  covGrad(1:size(covGrad, 1)+1:end) = 0;
%  base = covGradDiag.*diag(k)./(l.^3).*diag(n2);
%  for i = 1:size(g, 2) - 1
%    g(i) = g(i) + sum(base.*gOut(:, i));
%  end
%  base = covGrad.*k;
%  base = base.*(d/2*(L2.^4 - L1.^4) + 2*L1.*L1.*n2)./(w2.*w2.*L1);
%  for i = 1:size(g, 2)-1
%    g(i) = g(i) + sum(sum(base.*repmat(gOut(:, i), 1, size(x, 1))));
%  end
%~/
else
  [k, sk, n2, w2, l, l2] = gibbsKernCompute(kern, x, varargin{1});
  gOut = modelOutputGrad(kern.lengthScaleFunc, x);
  gradFact = fhandle(l, 'gradfact');
  gOut2 = modelOutputGrad(kern.lengthScaleFunc, varargin{1});
  gradFact2 = fhandle(l2, 'gradfact');
  L1 = repmat(l, 1, size(l2, 1));
  L2 = repmat(l2, 1, size(l, 1))';
  base = varargin{end}.*k;
  base2 = base.*(kern.inputDimension/2*(L2.^4 - L1.^4) + 2*L1.*L1.*n2)./(w2.*w2.*L1);
  for i = 1:size(g, 2)-1
    g(i) = g(i) + sum(sum(base2.*repmat(gOut(:, i).*gradFact, 1, size(varargin{1}, ...
                                                      1))));
  end
  base2 = base.*(kern.inputDimension/2*(L1.^4 - L2.^4) + 2*L2.*L2.* ...
                 n2)./(w2.*w2.*L2);
  for i = 1:size(g, 2)-1
    g(i) = g(i) + sum(sum(base2.*repmat(gOut2(:, i)'.*gradFact2', size(x, ...
                                                      1), 1)));
  end
  
end
g(end) =  sum(sum(varargin{end}.*sk));
