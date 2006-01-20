function [gmu, gsigmavar] = gpPosteriorGradMeanVar(model, X);

% GPPOSTERIORGRADMEANVAR Gadient of the mean and variances of the posterior at points given by X.

% FGPLVM


if ~isfield(model, 'alpha')
  model = gpComputeAlpha(model);
end

if size(X, 1) > 1
  error('This function only handles one data-point at a time')
end
switch model.approx
 case 'ftc'
  gX = kernGradX(model.kern, X, model.X);
  kX = kernCompute(model.kern, X, model.X)';
 case {'dtc', 'fitc', 'pitc'}
  gX = kernGradX(model.kern, X, model.X_u);
  kX = kernCompute(model.kern, X, model.X_u)';
 otherwise
  error('Unrecognised approximation type');
  
end
diaggK = kernDiagGradX(model.kern, X);

gmu = zeros(size(X, 2), model.d);
gsigmavar = zeros(size(X, 2), model.d);

switch model.approx
 case 'ftc'
  Kinvgk = model.invK_uu*gX;
 case 'dtc'
  Kinvgk = ((model.invK_uu - (1/model.beta)*model.Ainv)*gX);
 case {'fitc', 'pitc'}
  Kinvgk = (model.invK_uu - model.Ainv)*gX;
 otherwise
  error('Unrecognised approximation type');
end

gsigmavar = repmat(diaggK' - 2*Kinvgk'*kX, 1, model.d);
gmu = gX'*model.alpha; 

gmu = gmu.*repmat(model.scale, model.q, 1);
gsigmavar = gsigmavar.*repmat(model.scale, model.q, 1);
