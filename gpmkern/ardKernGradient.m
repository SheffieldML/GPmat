function g = ardKernGradient(kern, x, covGrad)

% ARDKERNGRADIENT Gradient of ARD kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% pre-built RBF and linear ARD
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
% SEEALSO : ardKernParamInit, kernGradient, ardKernDiagGradient, kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2004

% KERN


[k, rbfPart, linearPart, dist2xx] = ardKernCompute(kern, x);
g(1) = - .5*sum(sum(covGrad.*rbfPart.*dist2xx));
g(2) =  sum(sum(covGrad.*rbfPart))/kern.rbfVariance;
g(3) =  sum(sum(covGrad));
g(4) =  trace(covGrad);
g(5) =  sum(sum(covGrad.*linearPart))/kern.linearVariance;
for i = 1:size(x, 2)
  g(5+i) =  sum(sum(covGrad.*((kern.linearVariance)*x(:, i)*x(:, i)' ...
                                      -.5*(kern.inverseWidth)*dist2(x(:, i), ...
                                                    x(:, i)) ...
                                      .*rbfPart)));
end
