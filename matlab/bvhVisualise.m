function handle = bvhVisualise(channels, bvhStruct, padding)

% BVHVISUALISE For updating a bvh representation of 3-D data.

% MOCAP

if nargin<3
  padding = 0;
end
channels = [channels zeros(1, padding)];
vals = bvh2xyz(bvhStruct, channels);
connect = bvhConnectionMatrix(bvhStruct);

indices = find(connect);
[I, J] = ind2sub(size(connect), indices);
handle(1) = plot3(vals(:, 1), vals(:, 3), vals(:, 2), '.');
set(handle(1), 'markersize', 20);
%/~
%set(handle(1), 'visible', 'off')
%~/
hold on
grid on
for i = 1:length(indices)
  handle(i+1) = line([vals(I(i), 1) vals(J(i), 1)], ...
              [vals(I(i), 3) vals(J(i), 3)], ...
              [vals(I(i), 2) vals(J(i), 2)]);
  set(handle(i+1), 'linewidth', 2);
end
axis equal
xlabel('x')
ylabel('z')
zlabel('y')
axis on