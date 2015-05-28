function options = mvuOptions(neighbours)

% MVUOPTIONS Options for a MVU.
%
%	Description:
%
%	OPTIONS = MVUOPTIONS(NEIGHBOURS) returns the default options for
%	maximum variance unfolding.
%	 Returns:
%	  OPTIONS - default options structure for MVU.
%	 Arguments:
%	  NEIGHBOURS - the number of neighbours to use.
%	
%
%	See also
%	MVUCREATE, MODELCREATE


%	Copyright (c) 2009 Neil D. Lawrence


  if nargin < 1
    neighbours = 7;
  end
  options.numNeighbours = neighbours;
  options.solver = 1;
end