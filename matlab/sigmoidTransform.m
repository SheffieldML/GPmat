function y = sigmoidTransform(x, transform)

% SIGMOIDTRANSFORM Constrains a parameter to be between 0 and 1.

% IVM

limVal = 36;
y = zeros(size(x));
switch transform
 case 'atox'
  index = find(x<-limVal);
  y(index) = eps;
  x(index) = NaN;
  index = find(x<limVal);
  y(index) = sigmoid(x(index));
  x(index) = NaN;
  index = find(~isnan(x));
  y(index) = 1-eps;

 case 'xtoa'
  y = invSigmoid(x);
 case 'gradfact'
  y = x.*(1-x);
end
