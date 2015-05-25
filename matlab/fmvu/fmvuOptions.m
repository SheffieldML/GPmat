function options = fmvuOptions(neighbours)

% FMVUOPTIONS Create a default options structure for the FMVU model.
% FORMAT
% DESC creates a default options structure for the fast maximum variance unfolding model
% structure.
% ARG neighbours : the number of neighbours to use.
% RETURN options : the default options structure.
%
% SEEALSO : fmvuCreate, modelOptions
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  if nargin < 1
    neighbours = 7;
  end
  options.numNeighbours = neighbours;
end