function handle = xyzpoppeVisualise(pos,fid)

% XYZPOPPEVISUALISE Draw the Poppe figure return the graphics handle.
%
%	Description:
%
%	HANDLE = XYZPOPPEVISUALISE(POS, V) draws the stick figure from the
%	Poppe silhouette data.
%	 Returns:
%	  HANDLE - the graphics handle of the drawn object.
%	 Arguments:
%	  POS - the positions of the joints.
%	  V - the view point for the figure (defaults to standard 3D view).
%	
%
%	See also
%	XYZPOPPEDRAW, XYZPOPPEMODIFY


%	Copyright (c) 2008 Carl Henrik Ek and Neil Lawrence


if(nargin<2);
  fid = 1;
end

% Convert positions for plotting.
joint = xyzpoppe2joint(pos);
handle = xyzpoppeDraw(joint);
axis equal
view(3);

xlabel('x');
ylabel('z');
zlabel('y');
axis on

return;