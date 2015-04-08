function [y, a] = ivmfwd(x, model)

% IVMFWD Make out put predictions for the IVM

% IVM
D = size(model.y, 2);
kx = computeKernel(model.X(model.I, :), model.kern.lntheta, model.kern.type, x);
beta = zeros(length(model.I), length(model.Sigma));
for i = 1:D
  beta(:, i) = model.Sigma(i).Linv'*model.Sigma(i).Linv*model.m(model.I, i);
end

a = kx'*beta;
y = sign(a);

