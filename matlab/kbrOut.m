function Y = kbrOut(model, X);

% KBROUT Obtain the output of the kernel based regression model.

numData = size(X, 1);
Y = kernCompute(model.kern, X, model.X)*model.A+ones(numData, 1)*model.b;
