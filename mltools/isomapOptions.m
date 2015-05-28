function options = isomapOptions(neighbours)

% ISOMAPOPTIONS Options for a isomap.
% FORMAT
% DESC returns the default options for isomap.
% ARG neighbours : the number of neighbours to use.
% RETURN options : default options structure for isomap.
%
% SEEALSO : isomapCreate, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  if nargin < 1
    neighbours = 7;
  end
  options.numNeighbours = neighbours;
end
