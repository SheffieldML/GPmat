function [K1, K2] = kernelDiffParams(arg1, arg2, arg3);

% KERNELDIFFPARAMS Get gradients of kernel wrt its parameters.

if nargin < 3
  X1 = arg1;
  X2 = arg1;
  theta = arg2;
else
  X1 = arg1;
  X2 = arg2;
  theta = arg3;
end

theta = thetaConstrain(theta);

K = kernel(X1, X2, theta);
K2 = K/theta(2);
K1 = -0.5*dist2(X1, X2).*K;
