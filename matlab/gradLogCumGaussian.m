function y = gradLogCumGaussian(x)

% GRADLOGCUMGAUSSIAN Gradient of the log of the cumulative Gaussian.

% NDLUTIL

% Theoretically this is simply ngaussian(x)/cumGaussian(x) but there are
% problems in the tails of the distribution.

y = zeros(size(x));
index = find(x>0);
y(index) = ngaussian(x(index))./cumGaussian(x(index));
index = find(x<=0);
y(index) = 1./(sqrt(2*pi)*0.5*erfcx(-sqrt(2)/2*x(index)));
