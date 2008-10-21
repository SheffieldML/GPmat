function joint = xyzankur2joint(pos)

% XYZANKUR2JOINT
%
% COPYRIGHT : Carl Henrik Ek and Neil Lawrence, 2008

% MOCAP


joint(:,1) = pos(3:3:end);
joint(:,2) = -pos(1:3:end);
joint(:,3) = -pos(2:3:end);

return