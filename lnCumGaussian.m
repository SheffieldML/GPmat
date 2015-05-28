function y = lnCumGaussian(x)

% LNCUMGAUSSIAN log cumulative distribution for the normalised Gaussian.
% FORMAT
% DESC computes the logarithm of the cumulative Gaussian
% distribution.
% ARG X : input position.
% RETURN y : log probability of the value under the cumulative
% Gaussian.
%
% SEEALSO : erf, erfcx, cumGaussian, lnDiffCumGaussian, gaussOverDiffCumGaussian
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% NDLUTIL

index = find(x< 0);
if length(index)
  y(index) = -.5*x(index).*x(index) + log(.5) + log(erfcx(-sqrt(2)/2* ...
                                                    x(index)));
end
index = find(x>=0);
if length(index)
  y(index) = log(cumGaussian(x(index)));
end
y=reshape(y, size(x));
