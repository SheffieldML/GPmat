function options = leOptions(neighbours)

% LEOPTIONS Options for a Laplacian eigenmaps.
%
%	Description:
%
%	OPTIONS = LEOPTIONS(NEIGHBOURS) returns the default options for a
%	Laplacian eigenmaps.
%	 Returns:
%	  OPTIONS - default options structure for Laplacian eigenmaps.
%	 Arguments:
%	  NEIGHBOURS - the number of neighbours to use.
%	
%
%	See also
%	LECREATE, MODELCREATE


%	Copyright (c) 2009 Neil D. Lawrence


  if nargin < 1
    neighbours = 7;
  end
  options.numNeighbours = neighbours;
  options.isNormalised = true;
  options.weightType = 'rbf';
  options.weightScale = 1.0;
  options.regulariser = 0.0;
end
