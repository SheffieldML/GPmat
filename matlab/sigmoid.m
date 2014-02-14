function y = sigmoid(x)

% SIGMOID The sigmoid function

% GPMAT

y = ones(size(x))./(1+exp(-x));