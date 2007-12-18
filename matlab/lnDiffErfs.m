function v = lnDiffErfs(x1, x2),

% LNDIFFERFS Helper function for computing the log of difference
%   of two erfs.
% FORMAT
% DESC computes the log of the difference of two erfs in a numerically stable manner.
% ARG x1 : argument of the positive erf
% ARG x2 : argument of the negative erf
% RETURN v : log(erf(x1) - erf(x2))
%
% COPYRIGHT : Antti Honkela, 2007
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

% Case 1: arguments of different signs, no problems with loss of accuracy
I1 = sign(x1) ~= sign(x2);
% Case 2: both arguments are positive
% 2a: x1 > x2
I2a = (x1 > 0) & (x1 > x2) & ~I1;
% 2b: x1 <= x2
I2b = (x1 > 0) & ~I1 & ~I2a;
% Case 3: both arguments are negative
% 3a: x1 < x2
I3a = ~I1 & ~I2a & ~I2b & (x1 < x2);
% 3b: x1 >= x2
I3b = ~I1 & ~I2a & ~I2b & ~I3a;

warnState = warning('query', 'MATLAB:log:logOfZero');
warning('off', 'MATLAB:log:logOfZero');
v(I1) = log( erf(x1(I1)) - erf(x2(I1)) );
v(I2a) = log(erfcx(  x2(I2a)) ...
	     - erfcx(x1(I2a)) .* exp(x2(I2a).^2 - x1(I2a).^2)) ...
	 - x2(I2a).^2;
v(I2b) = log(erfcx(  x2(I2b)) .* exp(x1(I2b).^2 - x2(I2b).^2) ...
	     - erfcx(x1(I2b))) ...
	 - x1(I2b).^2;
v(I3a) = log(erfcx(  -x1(I3a)) .* exp(x2(I3a).^2 - x1(I3a).^2) ...
	     - erfcx(-x2(I3a))) ...
	 - x2(I3a).^2;
v(I3b) = log(erfcx(  -x1(I3b)) ...
	     - erfcx(-x2(I3b)) .* exp(x1(I3b).^2 - x2(I3b).^2)) ...
	 - x1(I3b).^2;
warning(warnState.state, 'MATLAB:log:logOfZero');
