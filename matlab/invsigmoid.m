function y=invsigmoid(x)

% INVSIGMOID The inverse of the sigmoid function.

% IVM

y = log(x./(1-x));