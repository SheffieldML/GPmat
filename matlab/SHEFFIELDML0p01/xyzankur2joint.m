function joint = xyzankur2joint(pos)

% XYZANKUR2JOINT Converts data to xyz positions for each joint.
%
%	Description:
%
%	XYZANKUR2JOINT(POS, JOINT) takes in a vector of values and returns a
%	matrix with points in rows and coordinate positions in columns.
%	 Arguments:
%	  POS - the vector of values.
%	  JOINT - the matrix of values with points in rows and x,y,z
%	   positions in columns.
%	


%	Copyright (c) 2008 Carl Henrik Ek and Neil Lawrence



joint(:,1) = pos(3:3:end);
joint(:,2) = -pos(1:3:end);
joint(:,3) = -pos(2:3:end);

return