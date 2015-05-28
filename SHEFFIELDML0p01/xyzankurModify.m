function handle = xyzankurModify(handle,pos, varargin)

% XYZANKURMODIFY  Helper function for modifying the point cloud from Agarwal and Triggs data.
%
%	Description:
%
%	HANDLE = XYZANKURMODIFY(HANDLE, JOINT) modifies the stick figure
%	associated with a set of x,y,z points for the Agarwal and Triggs
%	silhouette data.
%	 Returns:
%	  HANDLE - modified graphics handles.
%	 Arguments:
%	  HANDLE - graphics handles to modify.
%	  JOINT - the positions of the joints.
%	
%
%	See also
%	XYZANKURDRAW, XYZANKURVISUALISE


%	Copyright (c) 2008 Carl Henrik Ek and Neil Lawrence



joint = xyzankur2joint(pos);

if(iscell(handle))
  for(i = 1:1:length(handle))
    xyzankurDraw(joint,handle{i}); 
  end
else
  xyzankurDraw(joint,handle);
end
