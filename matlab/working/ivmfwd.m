function [y, a, p] = ivmfwd(x, model)

% IVMFWD Make out put predictions for the IVM

% IVM

kx = kernel(model.X(model.activeIndex, :), model.lntheta, model.kernelType, x);
Linv = eye(length(model.activeIndex))/model.L;
sqrtPi = diag(sqrt(model.sitePrecision(model.activeIndex))); 
betaTemp = sqrtPi*Linv';
beta = betaTemp*betaTemp'*model.siteMean(model.activeIndex);


a = kx'*beta;
y = sign(a);

