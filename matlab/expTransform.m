function y = expTransform(x, transform)

% EXPTRANSFORM Constrains a parameter to be positive through exponentiation.

% OPTIMI


limVal = 36;
y = zeros(size(x));
switch transform
 case 'atox'
  index = find(x<-limVal);
  y(index) = eps;
  x(index) = NaN;
  index = find(x<limVal);
  y(index) = exp(x(index));
  x(index) = NaN;
  index = find(~isnan(x));
  if ~isempty(index)
    y(index) = exp(limVal);
  end
 case 'xtoa'
  y = log(x);
 case 'gradfact'
  y = x;
end
