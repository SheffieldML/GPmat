function y = invCumGaussian(x)

% INVCUMGAUSSIAN Computes inverse of the cumulative Gaussian.
% FORMAT
% DESC computes the inverse of the cumulative Gaussian.
% ARG x : value between 0 and 1 to map onto the real line.
% RETURN y : the inverse of the cumulative Gaussian.
%
% SEEALSO : cumGaussian, erfinv
% 
% COPYRIGHT : Neil D. Lawrence, 2005

% NDLUTIL

y = erfinv(x*2 - 1)*2/sqrt(2);
