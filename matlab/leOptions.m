function options = leOptions(neighbours)

% LEOPTIONS Options for a Laplacian eigenmaps.
% FORMAT
% DESC returns the default options for a Laplacian eigenmaps.
% ARG neighbours : the number of neighbours to use.
% RETURN options : default options structure for Laplacian eigenmaps.
%
% SEEALSO : leCreate, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2009
  
% MLTOOLS

  if nargin < 1
    neighbours = 7;
  end
  options.numNeighbours = neighbours;
  options.isNormalised = true;
  options.weightType = 'rbf';
  options.weightScale = 1.0;
  options.regulariser = 0.0;
end
