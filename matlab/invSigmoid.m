function y=invSigmoid(x)

% INVSIGMOID The inverse of the sigmoid function.

% NDLUTIL

y = log(x./(1-x));