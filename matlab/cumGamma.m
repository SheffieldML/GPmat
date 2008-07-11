function y = cumGamma(x, a, b)

% CUMGAMMA Cumulative distribution for gamma.
% FORMAT
% DESC computes the cumulative gamma distribution.
% ARG x : input value.
% RETURN p : output probability.
%
% SEEALSO :  gammainc, gamma
%
% COPYRIGHT : Neil D. Lawrence, 2008

% NDLUTIL

y = gammainc(x*b, a);
