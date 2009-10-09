function options = pmvuOptions(neighbours)

% PMVUOPTIONS Create a default options structure for the PMVU model.
% FORMAT
% DESC creates a default options structure for the probabilistic maximum variance unfolding model
% structure.
% ARG neighbours : the number of neighbours to use.
% RETURN options : the default options structure.
%
% SEEALSO : pmvuCreate, modelOptions
%
% COPYRIGHT : Neil D. Lawrence 2009

% MLTOOLS

  if nargin < 1
    neighbours = 7;
  end
  options.numNeighbours = neighbours;
  options.isNormalised = false;
  options.regulariser = 0.0;
end
