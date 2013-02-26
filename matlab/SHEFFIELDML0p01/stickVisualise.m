function handle = stickVisualise(vals, connect)

% STICKVISUALISE For drawing a stick representation of 3-D data.
%
%	Description:
%
%	HANDLE = STICKVISUALISE(VALS, CONNECT) draws a stick man
%	representation in a 3-D plot.
%	 Returns:
%	  HANDLE - a vector of handles to the plotted structure.
%	 Arguments:
%	  VALS - the x,y,z channels to update the skeleton with.
%	  CONNECT - the connectivity of the skeleton.
%	
%
%	See also
%	STICKMODIFY


%	Copyright (c) 2005, 2006 Neil D. Lawrence


vals = reshape(vals, size(vals, 2)/3, 3);

indices = find(connect);
[I, J] = ind2sub(size(connect), indices);
handle(1) = plot3(vals(:, 1), vals(:, 2), vals(:, 3), '.');
set(handle(1), 'markersize', 20);
%set(handle(1), 'visible', 'off')
hold on
grid on
for i = 1:length(indices)
  handle(i+1) = line([vals(I(i), 1) vals(J(i), 1)], ...
              [vals(I(i), 2) vals(J(i), 2)], ...
              [vals(I(i), 3) vals(J(i), 3)]);
  set(handle(i+1), 'linewidth', 2);
end
axis equal
%set(gca, 'zlim', [-2 2])
%set(gca, 'ylim', [-2 2])
%set(gca, 'xlim', [-2 2])
%set(gca, 'cameraposition', [15.3758 -29.5366 9.54836])