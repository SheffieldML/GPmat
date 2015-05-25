function [X, sigma2] = isomapEmbed(Y, dims)

% ISOMAPEMBED Embed data set with Isomap.

% MLTOOLS

% Note: isomap code uses the transpose of a design matrix.
if any(any(isnan(Y)))
  error('Cannot initialise gplvm using isomap when missing data is present.')
else
  D = L2_distance(Y', Y', 1);
  options.dims = 1:dims;
  neighbours = 7;
  [Xstruct, sigma2, E] = Isomap(D, 'k', neighbours, options);
  X = zeros(size(Y, 1), 2);
  if length(Xstruct.index) ~= size(Y, 1)
    % We don't really deal with this problem correctly here ...
    warning('Isomap graph is not fully connected');
  end
  X(Xstruct.index, :) = Xstruct.coords{dims}';
  % Rescale X so that variance is 1 and mean is zero.
  meanX = mean(X);
  X = X-ones(size(Y, 1), 1)*meanX;
  varX = var(X);
  X = X*diag(sqrt(1./varX));
end
