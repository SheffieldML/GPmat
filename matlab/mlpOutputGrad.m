function g = mlpOutputGrad(model, X)

% MLPOUTPUTGRAD Evaluate derivatives of mlp model outputs with respect to parameters.

g = mlpderiv(model, X);