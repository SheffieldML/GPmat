function g = rbfOutputGrad(model, X)

% RBFOUTPUTGRAD Evaluate derivatives of rbf model outputs with respect to parameters.
% MLTOOLS

g = rbfderiv(model, X);