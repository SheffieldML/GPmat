function y = gradLogCumGaussian(x)

% GRADLOGCUMGAUSSIAN Gradient of the log of the cumulative Gaussian.

% IVM

% Theoretically this is simply ngaussian(x)/cumGaussian(x) but there are
% problems in the tails of the distribution.

    %.5*erfcx(-sqrt(2)/2*x)=exp(.5*x*x)*cumGaussian(x) ...

y = zeros(size(x));
index = find(x>0);
y(index) = ngaussian(x(index))./cumGaussian(x(index));
x(index) = NaN;
index = find(x<=0);
y(index) = 1./(sqrt(2*pi)*0.5*erfcx(-sqrt(2)/2*x(index)));
