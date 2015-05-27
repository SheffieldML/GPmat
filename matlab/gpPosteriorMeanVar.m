function [mu, varsigma] = gpPosteriorMeanVar(model, X);

% GPPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% FORMAT
% DESC returns the posterior mean and variance for a given set of
% points.
% ARG model : the model for which the posterior will be computed.
% ARG x : the input positions for which the posterior will be
% computed.
% RETURN mu : the mean of the posterior distribution.
% RETURN sigma : the variances of the posterior distributions.
%
% SEEALSO : gpCreate, gpPosteriorMeanCovar, gpPosteriorGradMeanVar
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006, 2009

% GP

model.m = gpComputeM(model);
if ~isfield(model, 'alpha')
  model = gpComputeAlpha(model);
end

maxMemory = 1000000;
switch model.approx
 case 'ftc'
  chunkSize = ceil(maxMemory/model.N);
 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
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
   case {'dtc', 'dtcvar', 'fitc', 'pitc'}
    KX_star = kernCompute(model.kern, model.X_u, X(indices, :));  
  end

  % Compute mean, using precomputed alpha vector.
  if ~isfield(model, 'isMissingData') || ~model.isMissingData ... 
    | ~strcmp(model.approx, 'ftc')
    mu(indices, :) = KX_star'*model.alpha;
  else
    for i = 1:model.d
      mu(indices, i) = KX_star(model.indexPresent{i}, :)' ...
          *model.alpha(model.indexPresent{i}, i);
    end
  end
  
  % Compute variances if requried.
  if nargout > 1
    if ~isfield(model, 'isSpherical') || model.isSpherical
      % Compute diagonal of kernel for new point.
      diagK = kernDiagCompute(model.kern, X(indices, :));
      switch model.approx
       case 'ftc'
        Kinvk = model.invK_uu*KX_star;
       case {'dtc', 'dtcvar', 'fitc', 'pitc'}
        Kinvk = (model.invK_uu - (1/model.beta)*model.Ainv)*KX_star;
      end
      varsig = diagK - sum(KX_star.*Kinvk, 1)';
      if isfield(model, 'beta')
        varsig = varsig + (1/model.beta);
      end
      varsigma(indices, :) = repmat(varsig, 1, model.d);
    else
      diagK = kernDiagCompute(model.kern, X(indices, :));
      for i = 1:model.d
        ind = model.indexPresent{i};
        switch model.approx
         case 'ftc'
          Kinvk = model.invK_uu{i}*KX_star(ind, :);
         otherwise 
          error(['Non-spherical not yet implemented for any approximation ' ...
                 'other than ''ftc''']);
        end
        varsigma(indices, i) = diagK - sum(KX_star(ind, :).*Kinvk, 1)';
      end
    end
  end
  
  
  % Rescale the mean
  mu(indices,:) = mu(indices, :).*repmat(model.scale, length(indices), 1);
  % Add the bias back in.
  mu(indices, :) = mu(indices, :) + repmat(model.bias, length(indices), 1);
  % If the mean function is present, add it it.
  if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
    mu(indices, :) = mu(indices, :) + modelOut(model.meanFunction, ...
                                               X(indices, :));
  end
  % rescale the variances
  if nargout > 1
    varsigma(indices, :) = varsigma(indices, :).* ...
        repmat(model.scale.*model.scale, length(indices), 1);
  end

  % Prepare for the next chunk.
  startVal = endVal + 1;
  endVal = endVal + chunkSize;
  if endVal > size(X, 1)
    endVal = size(X, 1);
  end
end
