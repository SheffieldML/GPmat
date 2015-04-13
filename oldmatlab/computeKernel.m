function [K, invK] = computeKernel(X, theta);

% COMPUTEKERNEL Compute the kernel matrix for data X with parameters theta.

theta = thetaConstrain(theta);
K = kernel(X, X, theta);
K = K + eye(size(X, 1))*1/theta(end);
  
if nargout > 1
  invK = pdinv(K);
end
