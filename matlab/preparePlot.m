function prepareplot(limitVal, figNo)

% PREPAREPLOT Helper function for tidying up the plot before printing.

% NDLUTIL

% GENERAL

axis equal

set(gca, 'xlim', limitVal);
set(gca, 'ylim', limitVal);
xlim = get(gca, 'xlim');
ylim = get(gca, 'ylim');

xlim(1) = floor(xlim(1));
xlim(2) = ceil(xlim(2));
ylim(1) = floor(ylim(1));
ylim(2) = ceil(ylim(2));

set(gca, 'xlim', xlim);
set(gca, 'ylim', ylim);

dirs(1) = max([xlim(1) ylim(1)]);
dirs(2) = min([xlim(2) ylim(2)]);

line(dirs, dirs);
line(dirs, dirs+1);
line(dirs, dirs-1);

zeroaxes(gca, 0.02, 18, 'times')

