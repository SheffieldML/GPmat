function y = sigmoid(x)

% SIGMOID The sigmoid function

% SHEFFIELDML

y = ones(size(x))./(1+exp(-x));