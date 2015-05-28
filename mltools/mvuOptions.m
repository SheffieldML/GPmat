function options = mvuOptions(neighbours)

% MVUOPTIONS Options for a MVU.
% FORMAT
% DESC returns the default options for maximum variance unfolding.
% ARG neighbours : the number of neighbours to use.
% RETURN options : default options structure for MVU.
%
% SEEALSO : mvuCreate, modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  if nargin < 1
    neighbours = 7;
  end
  options.numNeighbours = neighbours;
  options.solver = 1;
end
