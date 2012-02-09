function handle = xyzVisualise(xyzChannels, skel)

% XYZVISUALISE For drawing an xyz representation of 3-D data.
% FORMAT
% DESC draws a skeleton representation in a 3-D plot.
% ARG xyzChannels : the 3-D coordinate channels to update the skeleton with.
% ARG skel : the skeleton structure.
% RETURN handle : a vector of handles to the plotted structure.
%
% SEEALSO : skelVisualise
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% MODIFICATIONS : Alfredo A. Kalaitzis, 2012
  
% MOCAP

connect = skelConnectionMatrix(skel);

indices = find(connect);
[I, J] = ind2sub(size(connect), indices);
handle(1) = plot3(xyzChannels(:, 1), xyzChannels(:, 3), xyzChannels(:, 2), '.');
axis ij % make sure the left is on the left.
set(handle(1), 'markersize', 20);
%/~
%set(handle(1), 'visible', 'off')
%~/
hold on
grid on
for i = 1:length(indices)
  handle(i+1) = line([xyzChannels(I(i), 1) xyzChannels(J(i), 1)], ...
              [xyzChannels(I(i), 3) xyzChannels(J(i), 3)], ...
              [xyzChannels(I(i), 2) xyzChannels(J(i), 2)]);
  set(handle(i+1), 'linewidth', 2);
end
axis equal
xlabel('x')
ylabel('z')
zlabel('y')
axis on