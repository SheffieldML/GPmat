function y = sigmoidBound(x)

% SIGMOIDBOUND Constrains a parameter to be between 0 and 1.

limValue = 36;
index = find(x<-limValue);
y(index) = eps;
x(index) = NaN;
index = find(x<limValue);
y(index) = sigmoid(x(index));
x(index) = NaN;
index = find(~isnan(x));
y(index) = 1-eps;

