function handle = plot3Modify(handle, values, Y)

% PLOT3MODIFY Helper code for visualisation of 3-d data.
%
%	Description:
%	handle = plot3Modify(handle, values, Y)
%

set(handle, 'XData', values(1), 'YData', values(2), 'ZData', values(3));
