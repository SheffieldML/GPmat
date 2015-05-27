function handle = plot3Visualise(vals, Y)

% PLOT3VISUALISE  Helper code for plotting a plot3 visualisation.

% MLTOOLS

if length(vals)>3
  error('plot3Visualise requires output dimension of 3');
end

plot3(Y(:, 1), Y(:, 2), Y(:, 3), 'x'); hold on
handle = plot3(vals(1), vals(2), vals(3), 'o');
hold off
