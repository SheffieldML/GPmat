function [mu, varsigma] = gpPosteriorMeanVar(model, X);

% GPPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.

% FGPLVM


if ~isfield(model, 'alpha')
  model = gpComputeAlpha(model);
end

maxMemory = 1000000;
switch model.approx
 case 'ftc'
  chunkSize = ceil(maxMemory/model.N);
 case {'dtc', 'fitc', 'pitc'}
  chunkSize = ceil(maxMemory/model.k);
end

mu = zeros(size(X, 1), model.d);
if nargout > 1
  varsigma = zeros(size(X, 1), model.d);
end

startVal = 1;
endVal = chunkSize;
if endVal>size(X, 1)
  endVal = size(X, 1);
end

while startVal <= size(X, 1)
  indices = startVal:endVal;

  % Compute kernel for new point.
  switch model.approx
   case 'ftc'
    KX_star = kernCompute(model.kern, model.X, X(indices, :));  
   case {'dtc', 'fitc', 'pitc'}
    KX_star = kernCompute(model.kern, model.X_u, X(indices, :));  
  end

  % Compute mean, using precomputed alpha vector.
  mu(indices, :) = KX_star'*model.alpha;

  
  % Compute variances if requried.
  if nargout > 1
    % Compute diagonal of kernel for new point.
    diagK = kernDiagCompute(model.kern, X(indices, :));
    switch model.approx
     case 'ftc'
      Kinvk = model.invK_uu*KX_star;
     case 'dtc'
      Kinvk = ((model.invK_uu - model.sigma2*model.Ainv)*KX_star);
     case {'fitc', 'pitc'}
      Kinvk = (model.invK_uu - model.Ainv)*KX_star;
    end
    varsig = diagK - sum(KX_star.*Kinvk, 1)';
    if isfield(model, 'sigma2')
      varsig = varsig + model.sigma2;
    end
    varsigma(indices, :) = repmat(varsig, 1, model.d);
  end
  
  
  
  % Add the mean back in.
  mu(indices,:) = mu(indices, :).*repmat(model.scales, length(indices), 1) + repmat(model.m, length(indices), 1);
  
  % rescale the variances
  if nargout > 1
    varsigma(indices, :) = varsigma(indices, :).* ...
        repmat(model.scales.*model.scales, length(indices), 1);
  end

  % Prepare for the next chunk.
  startVal = endVal + 1;
  endVal = endVal + chunkSize;
  if endVal > size(X, 1)
    endVal = size(X, 1);
  end
end