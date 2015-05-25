function y = gradFactLinearBound(x)

% GRADFACTLINEARBOUND Gradient multiplier for linear bound.


limVal = 36;
index = find(x>limVal);
y(index) = 1;
index = find(x<=limVal);
y(index) = (exp(x(index))-1)./exp(x(index));
