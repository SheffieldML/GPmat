function [theta, rval] = do_distfit(data, dist, verbose),

p = [.05, .25, .50, .75, .95];

if nargin < 2,
  dist = @gaminv;
end

[theta, rval] = fminsearch(@(theta) distfit(data, theta, dist), [1, 1]);

if nargin > 2 && verbose,
  disp(data);
  disp(dist(p, theta(1), theta(2)));
end
