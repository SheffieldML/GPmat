function [gmu, gsigmavar, factors] = gpPosteriorGradMeanCovar(model, X);

% GPPOSTERIORGRADMEANCOVAR Gadient of the mean and variances of the posterior at points given by X.
% FORMAT
% DESC computes the gradient of the mean and covariances of the
% posterior distribution of a Gaussian process with respect to the
% input locations. 
% ARG model : the model for which gradients are to be computed.
% ARG X : the input locations where gradients are to be computed.
% RETURN gmu : the gradient of the posterior mean with respect to
% the input locations.
% RETURN gCovar : the gradients of the posterior covariance with
% respect to the input locations. By raw, we mean that the gradient
% has not yet been multiplied by any output scale in that direction
% (as is done for gpPosteriorGradMeanCovar). The gradients are
% stored in a cell array of dimension MODEL.q x MODEL.d. 
%
% DESC computes the gradient of the mean and covariances of the
% posterior distribution of a Gaussian process with respect to the
% input locations. Returns a compact representation for the
% covariances which separates the factors associated with the
% different dimensions from the covariance gradients.
% ARG model : the model for which gradients are to be computed.
% ARG X : the input locations where gradients are to be computed.
% RETURN gmu : the gradient of the posterior mean with respect to
% the input locations.
% RETURN grCovar : the 'raw' gradient of the posterior covariance with
% respect to the input locations. By raw, we mean that the gradient
% has not yet been multiplied by any output scale in that direction
% (as is done for gpPosteriorGradMeanCovar). The gradients are
% stored in a cell array of length MODEL.q. 
% RETURN factors : the factors for multiplying the 'raw' gradients
% of the covariances by. 
%
% SEEALSO : gpCreate, gpPosteriorMeanVar
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2009

% GP


if ~isfield(model, 'alpha')
  model = gpComputeAlpha(model);
end

switch model.approx
 case 'ftc'
  gX = kernGradX(model.kern, X, model.X);
  kX_star = kernCompute(model.kern, X, model.X)';
 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
  gX = kernGradX(model.kern, X, model.X_u);
  kX_star = kernCompute(model.kern, X, model.X_u)';
 otherwise
  error('Unrecognised approximation type');
  
end
K = kernGradX(model.kern, X);


if ~model.isMissingData
  for i = 1:model.q
    switch model.approx
     case 'ftc'
      KinvgK = model.invK_uu*squeeze(gX(:, i, :));
     case {'dtc', 'dtcvar', 'fitc', 'pitc'}
      KinvgK = (model.invK_uu - (1/model.beta)*model.Ainv)*squeeze(gX(:, i, :));
     otherwise
      error('Unrecognised approximation type');
    end
    kXTKinvgK = kX_star'*KinvgK;
    gCovar{i} = squeeze(K(:, i, :))-kXTKinvgK - diag(diag(kXTKinvgK));
    gmu{i} = squeeze(gX(:, i, :))'*model.alpha.*repmat(model.scale, ...
                                                      size(X, 1), 1);
  end
  
  
  % Deal with scaling.
  if nargout < 3
    for i = 1:model.q
      for j = 1:model.d
        gsigmavar{i, j} = gCovar{i}*model.scale(j)*model.scale(j);
      end
    end
  else
    factors = model.scale.*model.scale;
    gsigmavar = gCovar;
  end
else
  error('Not yet implemented for models trained on missing data.');
end
