function y = sigmoid(x)

% SIGMOID The sigmoid function

% IVM

y = ones(size(x))./(1+exp(-x));