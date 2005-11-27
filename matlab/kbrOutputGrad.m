function g = kbrOutputGrad(model, X)

% KBROUTPUTGRAD Evaluate derivatives of kernel based regression model outputs with respect to parameters.

numData = size(X, 1);
g(:, :, 1) = [model.K zeros(size(model.K)) ones(numData, 1) ...
                    zeros(numData, 1)];
g(:, :, 2) = [zeros(size(model.K)) model.K zeros(numData, 1) ...
                    ones(numData, 1)];
