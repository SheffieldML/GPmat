function y = negLogLogitTransform(x, transform)

% NEGLOGLOGIT Constrains a parameter to be positive.
%
% y = negLogLogitTransform(x, transform)
%

% Copyright (c) 2005 Neil D. Lawrence
% negLogLogitTransform.m version 1.1




y = zeros(size(x));
limVal = 36;
switch transform
 case 'atox'
  index = find(x<-limVal);
  y(index) = eps;
  x(index) = NaN;
  index = find(x<limVal);
  y(index) = log(1+exp(x(index)));
  x(index) = NaN;
  index = find(~isnan(x));
  y(index) = x(index);
 case 'xtoa'
  index = find(x<limVal);
  y(index) = log(exp(x(index))-1);
  index = find(x>=limVal);
  y(index) = x(index);
 case 'gradfact'
  index = find(x>limVal);
  y(index) = 1;
  index = find(x<=limVal);
  y(index) = (exp(x(index))-1)./exp(x(index));
end
  