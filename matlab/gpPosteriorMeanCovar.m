function [mu, covarsig, scales] = gpPosteriorMeanCovar(model, X);

% GPPOSTERIORMEANCOVAR Mean and covariances of the posterior at points given by X.
% FORMAT
% DESC gives the posterior mean and covariance at the points given
% by X.
% ARG model : the model for which the posterior will be computed.
% ARG X : the latent positions where the mean and covariance will
% be computed.
% RETURN mu : the posterior mean.
% RETURN Sigma : the posterior covariance *without scaling*.
%
% SEEALSO : gpCreate, gpPosteriorMeanVar
%
% COPYRIGHT : Neil D. Lawrence, 2006

% FGPLVM

mu = gpPosteriorMeanVar(model, X);

if nargout > 1
  if size(X, 1)>1000
    warning(['Computation of covariances takes a long time for larger ' ...
             'data sets, are you sure you did''nt just want ' ...
             'variances? If so use gpPosteriorMeanVar.'])
  end
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
    Kinvk = (model.invK_uu - (1/model.beta)*model.Ainv)*KX_star;
   case {'fitc', 'pitc'}
    Kinvk = (model.invK_uu - (1/model.beta)*model.Ainv)*KX_star;
  end
  
  covarsig = K - KX_star'*Kinvk;
  if isfield(model, 'beta')
    covarsig = covarsig + eye(size(X, 1))*(1/model.beta);
  end
end
    
