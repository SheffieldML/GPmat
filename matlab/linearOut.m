function Y = linearOut(model, X);

% LINEAROUT Obtain the output of the linear model.

numData = size(X, 1);
Y = X*model.W + ones(numData, 1)*model.b;
