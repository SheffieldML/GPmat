function bvhModify(handle, channels, skel, padding)

% BVHMODIFY Helper code for visualisation of bvh data.
%
%	Description:
%
%	BVHMODIFY(HANDLE, CHANNELS, SKEL) updates a bvh representation in a
%	3-D plot.
%	 Arguments:
%	  HANDLE - a vector of handles to the structure to be updated.
%	  CHANNELS - the channels to update the skeleton with.
%	  SKEL - the skeleton structure.
%	
%
%	See also
%	SKELMODIFY, BVHVISUALISE


%	Copyright (c) 2005, 2006 Neil D. Lawrence



if nargin<4
  padding = 0;
end

skelModify(handle, channels, skel, padding);