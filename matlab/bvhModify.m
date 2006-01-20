function bvhModify(handle, channels, skel, padding)

% BVHMODIFY Helper code for visualisation of bvh data.

% MOCAP


if nargin<4
  padding = 0;
end

skelModify(handle, channels, skel, padding);