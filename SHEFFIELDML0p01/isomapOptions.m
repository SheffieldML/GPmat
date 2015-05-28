function options = isomapOptions(neighbours)

% ISOMAPOPTIONS Options for a isomap.
%
%	Description:
%
%	OPTIONS = ISOMAPOPTIONS(NEIGHBOURS) returns the default options for
%	isomap.
%	 Returns:
%	  OPTIONS - default options structure for isomap.
%	 Arguments:
%	  NEIGHBOURS - the number of neighbours to use.
%	
%
%	See also
%	ISOMAPCREATE, MODELCREATE


%	Copyright (c) 2009 Neil D. Lawrence


  if nargin < 1
    neighbours = 7;
  end
  options.numNeighbours = neighbours;
end