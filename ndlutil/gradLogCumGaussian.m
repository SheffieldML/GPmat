function y = gradLogCumGaussian(x)

% GRADLOGCUMGAUSSIAN Gradient of the log of the cumulative Gaussian.
% FORMAT 
% DESC returns the gradient of the logarithm of the cumulative
% Gaussian. Theoretically this is simply ngaussian(x)/cumGaussian(x)
% but there are problems in the tails of the distribution. This
% function attempts to deal with these problems without numerical
% error creeping in.
% ARG x : the input to the function.
% ARG y : the gradient of the log cumulative Gaussian.
%
% SEEALSO : ngaussian, cumGaussian, erfcx
% 
% COPYRIGHT : Neil D. Lawrence, 2005

% NDLUTIL


y = zeros(size(x));
index = find(x>0);
y(index) = ngaussian(x(index))./cumGaussian(x(index));
index = find(x<=0);
y(index) = 1./(sqrt(2*pi)*0.5*erfcx(-sqrt(2)/2*x(index)));
