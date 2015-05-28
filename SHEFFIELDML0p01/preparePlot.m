function preparePlot(limitVal, ax)

% PREPAREPLOT Helper function for tidying up the plot before printing.
%
%	Description:
%
%	PREPAREPLOT(LIMITVAL, AX) is a helper function for tidying up a plot
%	before printing.
%	 Arguments:
%	  LIMITVAL - the limits to be applied to the axes.
%	  AX - the axes to apply the plot preparation to.
%	
%
%	See also
%	ZEROAXES


%	Copyright (c) 2005 Neil D. Lawrence


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

