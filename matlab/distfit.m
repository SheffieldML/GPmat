function [theta, rval] = distfit(data, dist, verbose),
% DISTFIT Fit a distribution to given parameter percentiles
%
%   [theta, rval] = DISTFIT(prctiles[, dist, verbose])
% returns the parameter vector theta and residual rval when fitting
% a distribution to the 5th, 25th, 50th, 75th and 95th percentile
% given in prctiles.
% The argument dist specifies the distribution to use. Possible
% values are:
%  'normal'  Gaussian, parametrised with mean and std     (default)
%  'gamma'   Gamma, parametrised in standard MATLAB style
% If the optional third argument is nonzero, the original percentiles
% and those of the fit are printed out in the end.
%
% COPYRIGHT : Antti Honkela, 2007
  
% SHEFFIELDML


if size(data, 1) > size(data, 2),
  data = data';
end

p = [.05, .25, .50, .75, .95];

if nargin < 2,
  dist = 'normal';
end

switch dist,
 case 'gamma',
  pdf = @gaminv;
 case 'normal',
  pdf = @norminv;
 otherwise,
  error('Unknown distribution.')
end

[theta, rval] = fminsearch(@(theta) distfit_obj(data, theta, pdf), [1, 1]);

if nargin > 2 && verbose,
  disp(data);
  disp(pdf(p, theta(1), theta(2)));
end


function r = distfit_obj(y, theta, pdf),

p = [.05, .25, .50, .75, .95];

x = pdf(p, theta(1), theta(2));
r = .5 * sum((x - y).^2);
