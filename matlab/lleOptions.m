function options = lleOptions(neighbours, latentDim)

% LLEOPTIONS Options for a density network.
% FORMAT
% DESC returns the default options for a locally linear embedding.
% RETURN options : default options structure for locally linear embedding.
%
% SEEALSO : lleCreate, mlpCreate, rbfCreate, kbrCreate
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

if nargin < 2
  latentDim = 2
  if nargin < 1
    neighbours = 7;
  end
end
options.latentDim = latentDim;
options.numNeighbours = neighbours;


