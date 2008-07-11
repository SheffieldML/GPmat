function preparePlot(limitVal, ax)

% PREPAREPLOT Helper function for tidying up the plot before printing.
% FORMAT
% DESC is a helper function for tidying up a plot before printing. 
% ARG limitVal : the limits to be applied to the axes.
% ARG ax : the axes to apply the plot preparation to.
% 
% COPYRIGHT : Neil D. Lawrence, 2005
% 
% SEEALSO : zeroAxes

% NDLUTIL

axis equal

if nargin < 2
  ax = gca;
end
set(ax, 'xlim', limitVal);
set(ax, 'ylim', limitVal);
xlim = get(ax, 'xlim');
ylim = get(ax, 'ylim');

xlim(1) = floor(xlim(1));
xlim(2) = ceil(xlim(2));
ylim(1) = floor(ylim(1));
ylim(2) = ceil(ylim(2));

set(ax, 'xlim', xlim);
set(ax, 'ylim', ylim);

dirs(1) = max([xlim(1) ylim(1)]);
dirs(2) = min([xlim(2) ylim(2)]);

line(dirs, dirs);
line(dirs, dirs+1);
line(dirs, dirs-1);

zeroaxes(ax, 0.02, 18, 'times')

