function options = lleOptions(neighbours)

% LLEOPTIONS Options for a locally linear embedding.
% FORMAT
% DESC returns the default options for a locally linear embedding.
% ARG neighbours : the number of neighbours to use.
% RETURN options : default options structure for locally linear embedding.
%
% SEEALSO : lleCreate, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2008, 2009

% MLTOOLS

  if nargin < 1
    neighbours = 7;
  end
  options.numNeighbours = neighbours;
  options.isNormalised = true;
  options.regulariser = 0.0;
  options.acyclic = false;
end
