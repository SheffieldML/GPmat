function y = invLinearBound(x)

% LINEARBOUND Constrains a parameter to be positive.

limVal = 36;
index = find(x<-limVal);
y(index) = eps;
x(index) = NaN;
index = find(x<limVal);
y(index) = log(1+exp(x(index)));
x(index) = NaN;
index = find(~isnan(x));
y(index) = x(index);
