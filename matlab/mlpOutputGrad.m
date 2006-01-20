function g = mlpOutputGrad(model, X)

% MLPOUTPUTGRAD Evaluate derivatives of mlp model outputs with respect to parameters.

% MLTOOLS

g = mlpderiv(model, X);