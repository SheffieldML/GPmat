function model = isomapOptimise(model, display, iters)

% ISOMAPOPTIMISE Optimise an ISOMAP model.
% FORMAT
% DESC optimises an isomap model.
% ARG model : the model to be optimised.
% RETURN model : the optimised model.
%
% SEEALSO : modelOptimise
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

  
  % Note: isomap code uses the transpose of a design matrix.
if any(any(isnan(model.Y)))
  error('Cannot run isomap when missing data is present.')
else
  D = L2_distance(model.Y', model.Y', 1);
  options.dims = 1:model.q;
  [Xstruct, sigma2, E] = Isomap(D, 'k', model.k, options);
  model.X = zeros(size(model.Y, 1), 2);
  if length(Xstruct.index) ~= size(model.Y, 1)
    % We don't really deal with this problem correctly here ...
    warning('Isomap graph is not fully connected');
  end
  model.X(Xstruct.index, :) = Xstruct.coords{model.q}';
  % Rescale X so that variance is 1 and mean is zero.
  meanX = mean(model.X);
  model.X = model.X-ones(model.N, 1)*meanX;
  varX = var(model.X);
  model.X = model.X*diag(sqrt(1./varX));
end
