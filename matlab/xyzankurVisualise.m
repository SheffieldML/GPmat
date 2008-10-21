function handle = xyzankurVisualise(pos,fid,v)

% XYZANKURVISUALISE
%
% COPYRIGHT : Carl Henrik Ek and Neil Lawrence, 2008

% MOCAP


if(nargin<2)
  fid = 1;
end

figure(fid);

joint = xyzankur2joint(pos);
handle = xyzankurDraw(joint);axis equal;view([1 0 0]);
set(gca,'XLim',[-20 20],'YLim',[-20 20],'ZLim',[-5 70]);

if(exist('v','var')&&length(v)==2)
  view(v(1),v(2));
else
  view([1 0 0])
end

return;