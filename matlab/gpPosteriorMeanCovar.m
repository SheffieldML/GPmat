function [mu, covarsigma] = gpPosteriorMeanCovar(model, X);

% GPPOSTERIORMEANCOVAR Mean and covariances of the posterior at points given by X.
%
% [mu, covarsigma] = gpPosteriorMeanCovar(model, X);
%

% Copyright (c) 2006 Neil D. Lawrence
% gpPosteriorMeanCovar.m version 1.2



mu = gpPosteriorMeanVar(model, X);

if nargout > 1
  covarsigma = zeros(size(X, 1), size(X, 1));
end
if size(X, 1)>1000
  warning(['Computation of covariances takes a long time for larger ' ...
           'data sets, are you sure you did''nt just want ' ...
           'variances?'])
end

% Compute kernel for new point.
switch model.approx
 case 'ftc'
  KX_star = kernCompute(model.kern, model.X, X);  
 case {'dtc', 'fitc', 'pitc'}
  KX_star = kernCompute(model.kern, model.X_u, X);  
end

% Compute covariances if requried.
if nargout > 1
  % Compute kernel for new point.
  K = kernCompute(model.kern, X);
  switch model.approx
   case 'ftc'
    Kinvk = model.invK_uu*KX_star;
   case 'dtc'
    Kinvk = ((model.invK_uu - (1/model.beta)*model.Ainv)*KX_star);
   case {'fitc', 'pitc'}
    Kinvk = (model.invK_uu - model.Ainv)*KX_star;
  end
  
  covarsig = K - KX_star'*Kinvk;
  if isfield(model, 'beta')
    covarsig = covarsig + eye(size(X, 1))*(1/model.beta);
  end
end


    
% rescale the variances
if nargout > 1
  if ~all(model.scale==1)
    for i = 1:model.d
      covarsigma{i} = covarsig*model.scale(i).*model.scale(i);
    end
  else 
    covarsigma = covarsig;
  end
end
