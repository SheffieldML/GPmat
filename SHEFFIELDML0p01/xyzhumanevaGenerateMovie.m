function xyzhumanevaGenerateMovie(Z,filename)

% XYZHUMANEVAGENERATEMOVIE
%
%	Description:
%	


%	Copyright (c) 2008 Carl Henrik Ek and Neil Lawrence



if(nargin<2)
  filename = 'humaneva_movie';
end

% front view
mov = avifile(strcat(filename,'_front.avi'));
h = xyzhumanevaVisualise(Z(1,:),1);
F = getframe(gca);
mov = addframe(mov,F);
for(i = 2:1:size(Z,1))
  xyzhumanevaModify(h,Z(i,:));
  F = getframe(gca);
  mov = addframe(mov,F);
end
mov = close(mov);
clf(gcf,'reset');
% isometric view
mov = avifile(strcat(filename,'_3d.avi'));
h = xyzhumanevaVisualise3d(Z(1,:),1);
F = getframe(gca);
mov = addframe(mov,F);
for(i = 2:1:size(Z,1))
  xyzhumanevaModify(h,Z(i,:));
  F = getframe(gca);
  mov = addframe(mov,F);
end
mov = close(mov);

return