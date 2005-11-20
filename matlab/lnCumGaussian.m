function y = lnCumGaussian(x)

% LNCUMGAUSSIAN log cumulative distribution for Gaussian.

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