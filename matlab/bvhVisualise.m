function handle = bvhVisualise(channels, skel, padding)

% BVHVISUALISE For updating a bvh representation of 3-D data.

% MOCAP

if nargin < 3
  padding = 0;
end
handle = skelVisualise(channels, skel, padding);