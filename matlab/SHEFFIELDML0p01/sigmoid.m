function y = sigmoid(x)

% SIGMOID The sigmoid function
%
%	Description:
%	y = sigmoid(x)
%

y = ones(size(x))./(1+exp(-x));