function y = sigmoid(x)

% SIGMOID The sigmoid function

% NDLUTIL

y = ones(size(x))./(1+exp(-x));
