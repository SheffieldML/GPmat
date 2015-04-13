function handle = stickVisualise(vals, connect)

% STICKVISUALISE For updateing a stick representation of 3-D data.

% MOCAP

vals = reshape(vals, size(vals, 2)/3, 3);

indices = find(connect);
[I, J] = ind2sub(size(connect), indices);
handle(1) = plot3(vals(:, 1), vals(:, 2), vals(:, 3), '.');
set(handle(1), 'markersize', 12);
set(handle(1), 'visible', 'off')
hold on
grid on
for i = 1:length(indices)
  handle(i+1) = line([vals(I(i), 1) vals(J(i), 1)], ...
              [vals(I(i), 2) vals(J(i), 2)], ...
              [vals(I(i), 3) vals(J(i), 3)]);
  set(handle(i+1), 'linewidth', 2);
end
axis equal
set(gca, 'zlim', [-2 2])
set(gca, 'ylim', [-2 2])
set(gca, 'xlim', [-2 2])
set(gca, 'cameraposition', [15.3758 -29.5366 9.54836])