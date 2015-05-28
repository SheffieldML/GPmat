function [v, signs] = lnDiffErfs(x1, x2),

% LNDIFFERFS Helper function for computing the log of difference
%   of two erfs.
% FORMAT
% DESC computes the log of the difference of two erfs in a numerically stable manner.
% ARG x1 : argument of the positive erf
% ARG x2 : argument of the negative erf
% RETURN v : log(abs(erf(x1) - erf(x2)))
% RETURN s : sign(erf(x1) - erf(x2))
%
% FORMAT
% DESC computes the log of the difference of two erfs in a numerically stable manner.
% ARG x1 : argument of the positive erf
% ARG x2 : argument of the negative erf
% RETURN v : log(erf(x1) - erf(x2))     (Can be complex)
%
% COPYRIGHT : Antti Honkela, 2007, 2008
%
% MODIFICATIONS : David Luengo, 2009
%
% SEEALSO : gradLnDiffErfs

% NDLUTIL

x1 = real(x1);
x2 = real(x2);

v = zeros(max(size(x1), size(x2)));

if numel(x1) == 1,
  x1 = x1 * ones(size(x2));
end

if numel(x2) == 1,
  x2 = x2 * ones(size(x1));
end

signs = sign(x1 - x2);
I = signs == -1;
swap = x1(I);
x1(I) = x2(I);
x2(I) = swap;

% Case 1: arguments of different signs, no problems with loss of accuracy
I1 = (x1.*x2)<0;
% Case 2: x1 = x2
I2 = x1 == x2;
% Case 3: Both arguments are non-negative
I3 = (x1 > 0) & ~I1 & ~I2;
% Case 4: Both arguments are non-positive
I4 = ~I1 & ~I2 & ~I3;

warnState = warning('query', 'MATLAB:log:logOfZero');
warning('off', 'MATLAB:log:logOfZero');
v(I1) = log( erf(x1(I1)) - erf(x2(I1)) );
v(I2) = -inf;
v(I3) = log(erfcx(  x2(I3)) ...
	    - erfcx(x1(I3)) .* exp(x2(I3).^2 - x1(I3).^2)) ...
	- x2(I3).^2;
v(I4) = log(erfcx(  -x1(I4)) ...
	    - erfcx(-x2(I4)) .* exp(x1(I4).^2 - x2(I4).^2)) ...
	- x1(I4).^2;
warning(warnState.state, 'MATLAB:log:logOfZero');

if nargout < 2,
  v(I) = v(I) + pi*1i;
end
