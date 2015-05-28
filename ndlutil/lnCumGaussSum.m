function y = lnCumGaussSum(u1, u2, w1, w2)

% LNCUMGAUSSSUM The log of the weighted sum of two cumulative Gaussians.
% FORMAT
% DESC returns the logarithm of the weighted sum of two cumulative
% Gaussians.
% ARG u1 : argument of the first cumulative Gaussian.
% ARG u2 : argument of the second cumulative Gaussian.
% ARG w1 : weight of the first cumulative Gaussian.
% ARG w2 : weight of the second cumulative Gaussian.
%
% SEEALSO : cumGaussian, lnCumGaussian, lnDiffCumGaussian
%
% COPYRIGHT : Neil D. Lawrence, 2004

% NDLUTIL

y = zeros(size(u1));
safeCond = u1 > 0 & u2 > 0;
index = find(safeCond);
if ~isempty(index)
  y(index) = log(w1.*cumGaussian(u1(index)) ...
                 + w2*cumGaussian(u2(index)));
end
index = find(~safeCond & u1>u2);
if ~isempty(index)
  y(index) = log(w1) + lnCumGaussian(u1(index))...
      + log(1 + w2/w1*exp(lnCumGaussian(u2(index))...
                          -lnCumGaussian(u1(index))));
end
index = find(~safeCond & u2>=u1);
if ~isempty(index)
  y(index) = log(w2) + lnCumGaussian(u2(index))...
      + log(1 + w1/w2*exp(lnCumGaussian(u1(index))...
                          -lnCumGaussian(u2(index))));
end
