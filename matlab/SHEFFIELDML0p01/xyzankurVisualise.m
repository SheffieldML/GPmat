function handle = xyzankurVisualise(pos, fid, v)

% XYZANKURVISUALISE Draw the Agarwal & Triggs figure return the graphics handle.
%
%	Description:
%
%	HANDLE = XYZANKURVISUALISE(POS, V) draws the stick figure from the
%	Agarwal and Triggs silhouette data.
%	 Returns:
%	  HANDLE - the graphics handle of the drawn object.
%	 Arguments:
%	  POS - the positions of the joints.
%	  V - the view point for the figure (defaults to standard 3D view).
%	
%	
%
%	See also
%	XYZANKURDRAW, XYZANKURMODIFY


%	Copyright (c) 2008 Carl Henrik Ek and Neil Lawrence



if(nargin<2)
  fid = 1;
end

figure(fid);

% Convert positions for plotting.
joint = xyzankur2joint(pos);
handle = xyzankurDraw(joint);
axis equal
set(gca,'XLim',[-15 15],'YLim',[-15 15],'ZLim',[0 70]);
view(3)

if(exist('v','var')&&length(v)==2)
  view(v(1),v(2));
end
xlabel('x')
ylabel('z')
zlabel('y')
axis on

return;