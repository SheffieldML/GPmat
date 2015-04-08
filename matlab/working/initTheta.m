function lntheta = initTheta(kernelType, numIn)

% INITTHETA Initialise the kernel parameters.

% IVM

switch kernelType
 case 'linear'
  lntheta = [zeros(1, 3)];
 case 'rbf'
  lntheta = [zeros(1, 4)];
 case 'regular'
  lntheta = [zeros(1, 5)];
 case 'ARD'
  lntheta = [zeros(1, 5) zeros(1, numIn)];
end
