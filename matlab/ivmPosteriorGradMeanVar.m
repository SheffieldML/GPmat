function [gmu, gsigmavar] = ivmPosteriorGradMeanVar(model, X);

% IVMPOSTERIORGRADMEANVAR Gradient of mean and variances of the posterior wrt X.

% IVM

D = size(model.y, 2);
if size(X, 1) > 1
  error('This function only handles one data-point at a time')
end

gX = kernGradX(model.kern, X, model.X(model.I, :));
kX = kernCompute(model.kern, X, model.X(model.I, :))';
diaggK = kernDiagGradX(model.kern, X);


gmu = zeros(size(X, 2), D);
gsigmavar = zeros(size(X, 2), D);

if length(model.Sigma)>1
  for i = 1:D
    Kinvgk = model.Sigma(i).Linv'*(model.Sigma(i).Linv*gX);
    gsigmavar(:, i) = diaggK' - 2*Kinvgk'*kX;
    gmu(:, i) = Kinvgk'*model.m(model.I, i); 
  end 
else
  Kinvgk = model.Sigma.Linv'*(model.Sigma.Linv*gX);
  gsigmavar = repmat(diaggK' - 2*Kinvgk'*kX, 1, D);
  gmu = Kinvgk'*full(model.m(model.I, :)); 
end
