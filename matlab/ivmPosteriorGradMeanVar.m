function [gmu, gsigmavar] = ivmPosteriorGradMeanVar(X, model);

% IVMPOSTERIORGRADMEANVAR Gradient of mean and variances of the posterior wrt X.

% IVM

D = size(model.y, 2);
if size(X, 1) > 1
  error('This function only handles one data-point at a time')
end
gX = computeGradKernel(X, ...
                   model.kern.lntheta, ...
                   model.kern.type, ...
                   model.X(model.I, :));
kX = computeKernel(X, ...
                   model.kern.lntheta, ...
                   model.kern.type, ...
                   model.X(model.I, :))';

diaggK = kernelGradDiag(X, model.kern.lntheta, model.kern.type);

gmu = zeros(size(X, 2), D);
gsigmavar = zeros(size(X, 2), D);
if strcmp(model.noise.type, 'gaussian')
%  Lk = model.Sigma.Linv*kX;
%  Kinvk = model.Sigma.Linv'*Lk;
  Kinvgk = model.Sigma.Linv'*model.Sigma.Linv*gX;
  gsigmavar = repmat(diaggK' - 2*Kinvgk'*kX, 1, D);
end
for i = 1:D
  if ~strcmp(model.noise.type, 'gaussian')
%    Lk = model.Sigma(i).Linv*kX;
%    Kinvk = model.Sigma(i).Linv'*Lk;
    Kinvgk = model.Sigma(i).Linv'*model.Sigma(i).Linv*gX;
    gsigmavar(:, i) = diaggK' - 2*Kinvgk'*kX;
  end 
  gmu(:, i) = Kinvgk'*model.m(model.I, i); 
end